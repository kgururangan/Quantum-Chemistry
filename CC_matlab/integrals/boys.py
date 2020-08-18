#!/usr/bin/env python3

from scipy.special import hyp1f1
import numpy as np
import argparse

def boys(n,T):
    return hyp1f1(n+0.5,n+1.5,-T)/(2.0*n+1.0)

def main(args):
    n = args.n
    T = args.T
    print('F_{}({}) = {}'.format(n,T,boys(n,T)))
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate Boys function F_n(T) using Scipy')
    #parser.add_argument('-n','--n',type=int,help='integer n of Boys function')
    #parser.add_argument('-n','--n',type=int,help='integer n of Boys function')
    parser.add_argument('n',type=int,help='order n of Boys function')
    parser.add_argument('T',type=float,help='argument T of Boys function')
    args = parser.parse_args()
    main(args)
