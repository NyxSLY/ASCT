import subprocess as sp
from os import listdir
import logging
import sys
import multiprocessing as mp
import re
import os


logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)


def redup_then2bed(path):
    cmds = []
    for f in listdir(path):
        if f.find(".bam") != -1:
            in_file = path + '/' + f
            out_file = path + '/redup/' + re.sub('.bam','.redup.bam',f)
            metrics_file = path + '/redup/' + re.sub('.bam*','.redup.metrics',f)
            bed_out = "/".join(path.split("/")[:-1]) + '/bed/' + str(f).split(".")[0] + '.redup.bed'
            cmd = "java -jar picard.jar MarkDuplicates INPUT={} OUTPUT={} METRICS_FILE={} REMOVE_DUPLICATES=true".format(
                in_file, out_file, metrics_file)
            #sp.call(cmd, shell=True)
            # convert to bed format
            cmd1 = 'bedtools bamtobed -i {} > {}'.format(out_file, bed_out)
            #sp.call(cmd1, shell=True)
            cmds.append((cmd, cmd1))
    return cmds


def run_cmds(cmds):  # cmds should be a turple
    for cmd in cmds:
        logging.info('Started: {}'.format(cmd))
        see = sp.check_call(cmd, shell=True)
        if see == 0:
            pass
        else:
            logging.warning('shell command:\n{} has error!!\nscript stopped!')
            sys.exit()


def multi_run(func, argvs):  # func is the function defined by myself. argvs is a turple of turples
    processes = [mp.Process(target=func, args=(x,)) for x in argvs]
    for p in processes:
        p.start()

    for p in processes:
        p.join()


if __name__ == '__main__':
    try:
        os.mkdir('redup')
    except:
        pass
    try:
        os.mkdir('bed')
    except:
        pass
    cmds = redup_then2bed(sys.argv[1])
    #multi_run(run_cmds, cmds)
    run_cmds(cmds)
