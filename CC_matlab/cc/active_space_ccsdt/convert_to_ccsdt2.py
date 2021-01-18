import argparse
import sys
import os

def is_included(string):
    outside_projs = ['T3A.PPPhhH','T3A.PPPhhh','T3A.PPphhH','T3A.PPphhh', \
                     'T3A.PppHHH','T3A.PpphHH','T3A.PpphhH','T3A.Ppphhh', \
                     'T3A.pppHHH','T3A.ppphHH','T3A.ppphhH','T3A.ppphhh', \
                     'T3B.PPPhHh','T3B.PPPhhH','T3B.PPPhhh', \
                     'T3B.PpPhHh','T3B.PpPhhH','T3B.PpPhhh', \
                     'T3B.PPphHh','T3B.PPphhH','T3B.PPphhh', \
                     'T3B.PppHHH','T3B.PpphHH','T3B.PppHHh', \
                     'T3B.PpphHh','T3B.PpphhH','T3B.Ppphhh', \
                     'T3B.ppPHHH','T3B.ppPhHH','T3B.ppPHHh', \
                     'T3B.ppPhHh','T3B.ppPhhH','T3B.ppPhhh', \
                     'T3B.pppHHH','T3B.ppphHH','T3B.pppHHh', \
                     'T3B.ppphHh','T3B.ppphhH','T3B.ppphhh', \
                     'T3C.PPPhhH','T3C.PPPHhh','T3C.PPPhhh', \
                     'T3C.pPPhhH','T3C.pPPHhh','T3C.pPPhhh', \
                     'T3C.PPphhH','T3C.PPpHhh','T3C.PPphhh', \
                     'T3C.pPpHHH','T3C.pPphHH','T3C.pPpHhH', \
                     'T3C.pPphhH','T3C.pPpHhh','T3C.pPphhh', \
                     'T3C.PppHHH','T3C.PpphHH','T3C.PppHhH', \
                     'T3C.PpphhH','T3C.PppHhh','T3C.Ppphhh', \
                     'T3C.pppHHH','T3C.ppphHH','T3C.pppHhH', \
                     'T3C.ppphhH','T3C.pppHhh','T3C.ppphhh', \
                     'T3D.PPPhhH','T3D.PPPhhh','T3D.PPphhH','T3D.PPphhh', \
                     'T3D.PppHHH','T3D.PpphHH','T3D.PpphhH','T3D.Ppphhh', \
                     'T3D.pppHHH','T3D.ppphHH','T3D.ppphhH','T3D.ppphhh']

    if 'T3' not in string:
        return True
    else:
        if any([x in string for x in outside_projs]):
            return False
        else:
            return True

def clean_line(line):
    L = line.strip().split()
    cleaned_line = [x.strip() for x in L if is_included(x)]
    #temp = cleaned_line.copy()
    #for element in temp:
    #    print(element)
    #    if not element:
    #        cleaned_line.remove(element)
    return ' '.join(cleaned_line)


def main(args):
    f_in = args.f_in
    fname = f_in.split('.')
    f_target = fname[0]+'_ccsdt2.m'
    f_out = open(f_target,'w')

    with open(f_in,'r') as f:
        for line in f.readlines():
            cleaned_line = clean_line(line)
            f_out.writelines(cleaned_line)
            f_out.write('\n')
            #for element in cleaned_line.split():
            #    print(element)
            #    if element:
            #        f_out.write(element)
            #    else:
            #        continue
            #f_out.write('\n')
    f_out.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parser to convert fully active-spaced decomposed updates into CCSDt(II) form by removing lines involving T3 amplitudes that do not belong to the relevant projection space.')

    parser.add_argument('f_in',type=str,help='Input .m file containing full update')

    args = parser.parse_args()
    main(args)
