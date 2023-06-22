#!/bin/bash

name=$1

g++ --std=c++14 "$1".cpp -o $1 -Wall -O3
