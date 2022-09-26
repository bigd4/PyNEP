"""
Shuffle XYZ
===========================

This example shows how to shuffle structures randomly

provided by Kick-H
"""
from pynep.io import load_nep,dump_nep
import random

train_ratio = 0.8
rand = True

train_data = load_nep('train.in')

nframes = [i for i in range(len(train_data))]
if rand: random.shuffle(nframes)

train_frame = nframes[:int(len(train_data)*train_ratio)]
dump_nep('train.xyz', [train_data[i] for i in train_frame], ftype="exyz")

test_frame = nframes[int(len(train_data)*train_ratio):]
dump_nep('test.xyz', [train_data[i] for i in test_frame], ftype="exyz")
