#!/usr/bin/env python3

import os
import os.path as osp

import numpy as np
import pandas as pd

import torch as t
from torch.autograd import Variable

import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, random_split, DataLoader

import torchvision
import torchvision.transforms as transforms

# TODO : Implement transform to make n_transform x n_origial_samples dataset

class SpotDataset(Dataset):
    def __init__(self,
                 data_pth,
                 transform = None,
                 genes = None,
                 size = 1000,
                 n_classes = 5,
                 encoding = None,
                 arraysize = 8,
                ):

        self.data_pth = data_pth
        self.genelist = pd.Index(genes) if genes else pd.Index([])
        self.genes = genes
        self.samples = []
        self.arraysize = arraysize
        self.size = size
        self.nclasses = n_classes

        self.encoding =  {'her2lum' : 0,
                          'her2nonlum' : 1,
                          'luma' : 2,
                          'lumb' : 3,
                          'tnbc' : 4
                         }

        self.transform = transform

        self._init()

    def __len__(self,):
        return len(self.samples)

    def __getitem__(self,idx):
        return self.samples[idx]



    def _init(self,):
        count_pth = osp.join(self.data_pth,
                             'count_data',
                            )
        label_pth = osp.join(self.data_pth,
                            'label_data'
                            )

        count_pth = [osp.join(count_pth,x) for x in os.listdir(count_pth)]
        label_pth = [osp.join(label_pth,x) for x in os.listdir(label_pth)]

        count_pth.sort()
        label_pth.sort()

        self.nsamples = np.min((len(count_pth),self.size))

        unsrt = np.random.choice(np.arange(0,len(count_pth)),
                                 nsamples,
                                 replace = False
                                )

        count_pth = [ count_pth[x] for x in unsrt ]
        label_pth = [ label_pth[x] for x in unsrt ]

        for pcnt, plbl in zip(count_pth,label_pth):
            c = pd.read_csv(pcnt,
                            sep = '\t',
                            header = 0,
                            index_col = 0
                           )

            if self.genes is not None:
                tc = pd.DataFrame(np.zeros((self.samples[ii][0].shape[0],
                                            len(self.genelist))
                                           ),
                                 columns = self.genelist
                                 )

                inter = c.columns.intersection(self.genelist)
                tc.loc[:,inter] = c.values.astype(np.float32)
                tc = self.frame2tensor(tc)
                #if self.transform:
                #    tc = self.transform(tc)
            else:
                self.genelist = self.genelist.union(c.columns)
                tc = c

            with open(plbl,'r+') as lopen:
                l = lopen.readlines()[0]

            l = self._one_hot(l)

            self.samples.append([tc,l])

        if self.genes is None:
            for ii in range(len(self.samples)):

                tc = pd.DataFrame(np.zeros((self.samples[ii][0].shape[0],
                                            len(self.genelist))
                                           ),
                                 columns = self.genelist
                                 )

                inter = self.samples[ii][0].columns.intersection(self.genelist)
                tc.loc[:,inter] = self.samples[ii][0].values.astype(np.float32)
                self.samples[ii][0] = self.frame2tensor(tc.values.astype(np.float32))
                #if self.transform:
                #    self.samples[ii][0] = self.transform(self.samples[ii][0])

        self.samples = [ tuple(x) for x in self.samples]

        self.G = len(self.genelist)


    def _one_hot(self,x : str):

        return self.encoding[x]

    def frame2tensor(self,x):
        if not isinstance(x,np.ndarray):
            x = np.asarray(x)

        x = x.reshape(len(self.genelist),
                      self.arraysize,
                      self.arraysize,
                      )

        return x


class CNN(t.nn.Module):

    def  __init__(self,
                  in_channels):

        super(CNN, self).__init__()

        self.in_channels = in_channels
        # Input 8x8xG --> Output 8x8x100
        self.conv1 = nn.Conv2d(in_channels = self.in_channels,
                               out_channels = 100,
                               kernel_size = 3,
                               padding = 1,
                               stride = 1,
                              )

        # Input 8x8x100 --> Output 4x4x100
        self.pool = nn.AvgPool2d(kernel_size = 2,
                                 padding = 0,
                                 stride = 2,
                                )


        # Input 4x4x100 --> Output 3x3x20

        self.conv2 = nn.Conv2d(in_channels = 100,
                               out_channels = 20,
                               kernel_size  = 2,
                               padding = 0,
                               stride = 1,
                              )
        # Input 180x1 -- > Output 64
        self.fc1 = nn.Linear(in_features = 180,
                             out_features = 64,
                             bias = True,
                            )

        self.fc2 = nn.Linear(in_features = 64,
                             out_features = 5,
                             bias = True,
                            )

    def forward(self,x):

        x = F.relu(self.conv1(x))
        x = self.pool(x)
        x = self.conv2(x)
        x = x.view(-1,180)
        x = F.relu(self.fc1(x))
        x = self.fc2(x)

        return x


def train(net,
          train_set,
          val_set,
          batch_size,
          n_epochs,
          lr = 0.001,
          n_prints = 10,
          ):

    train_loader = DataLoader(train_set,
                              batch_size = batch_size,
                              num_workers = 2
                             )


    n_batches = len(train_loader)

    loss_fun = nn.CrossEntropyLoss()

    optim = t.optim.Adam(net.parameters(),
                         lr = lr,
                        )

    print_on = 10
    net.train()

    for epoch in range(n_epochs):

        total_loss = 0.0

        for i, data in enumerate(train_loader):

            inputs, labels = data
            inputs, labels = Variable(inputs), Variable(labels)

            optim.zero_grad()

            outputs = net(inputs)
            loss = loss_fun(outputs, labels.type(t.LongTensor))
            loss.backward()
            optim.step()

            total_loss += loss.item()

            if (i + 1) % print_on == 0:
                print(f"Epoch : {epoch + 1:d} | train_loss : {total_loss / print_on}")


        total_val_loss = 0.0
        val_loader = DataLoader(val_set,
                                  batch_size = batch_size,
                                  num_workers = 2
                                  )

        net.eval()

        for inputs, labels in val_loader:

             inputs, labels = Variable(inputs), Variable(labels)
             val_outputs = net(inputs)
             val_loss = loss_fun(val_outputs, labels.type(t.LongTensor))
             total_val_loss += val_loss.item()

        print(f"Validation loss : {total_val_loss / len(val_loader) }")

    print("Training Completed")

degs =[0,90,180,270]
p_train = 0.8
nsamples = 400

data_pth = "/home/alma/ST-2018/BC_Stanford/data/fakedata"

trf = transforms.Compose([transforms.ToTensor(),
                          transforms.RandomVerticalFlip(),
                          transforms.RandomHorizontalFlip(),
#                          transforms.Normalize(),
#                          transforms.RandomRotation(degrees = degs), # implement own rotation trans
                         ]
                        )


dataset = SpotDataset(data_pth,
                      size = nsamples,
                      transform = trf)


len_train = int(len(dataset)*p_train)
len_val = len(dataset) - len_train

train_set, val_set = random_split(dataset,
                                  (len_train,len_val),
                                  )


cnn_net = CNN(dataset.G)

train(cnn_net,
      train_set = train_set,
      val_set = val_set,
      batch_size = 32,
      n_epochs = 100,
      lr = 0.001,
     )






