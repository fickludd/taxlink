#!/usr/bin/python

import sys
import os
from pandas import *
import scipy.stats


dryrun = False
dual = False
usage = "usage:\n> taxlink2PyProphet [-d] [--dual] config.conf target.csv decoy.csv"

def readArgs(args):
	if len(args) == 3:
		global confPath
		global targetPath
		global decoyPath
		confPath 	= args[-3]
		targetPath 	= args[-2]
		decoyPath 	= args[-1]
	elif args[0] == "-d":
		global dryrun
		dryrun = True
		readArgs(args[1:])
	elif args[0] == "--dual":
		global dual
		dual = True
		readArgs(args[1:])
		

if len(sys.argv) < 4:
	print usage
	exit(1)
else:	
	readArgs(sys.argv[1:])


# CONF READING
def readConf(path, dryrun = True):
	conf = {}
	f = open(path, "r")
	for line in f:
		if '=' in line:
			kv = line.split('=', 2)
			conf[kv[0].strip()] = kv[1].strip()
	f.close()
	
	if dryrun:
		for key in conf:
			print key, "->", conf[key]
	
	return conf


def require(conf, name):
	if name not in conf:
		print "Required parameter '%s' not found in conf file '%s'. Exiting" % (name, confPath)
		exit(2)
	else:
		return conf[name]

		
def default(conf, name, x):
	if name not in conf:
		return x
	else:
		return conf[name]

# FILE PATH MANIPULATION
def withoutExt(f):
	return ".".join(f.split(".")[:-1])

def base(f):
	return withoutExt(os.path.basename(f))

def to(f, ext):
	return "%s.%s" % (withoutExt(f), ext)
	

conf = readConf(confPath)


if dryrun:
	print "  realFile:", targetPath
	print " decoyFile:", decoyPath
	print "      conf:", confPath
	for k,v in conf.iteritems():
		print "          ", k, ":", v
	
	

def setup(path, isDecoy):
	df = read_csv(path, sep="\t")
	df['transition_group_id'] = df['id']
	df['main_var_light_rt_score'] = df['light_rt_score']
	if dual:
		df['var_heavy_rt_score'] = df['heavy_rt_score']
		df['var_double_int_score'] = df['double_int_score']
		nXL = 2
	else:
		df['var_single_int_score'] = df['single_int_score']
		nXL = 1
	df['var_rel_mz_diff'] = abs(df['target_mz'] - df['mz']) * 2.0 / (df['target_mz'] + df['mz'])
	df['var_openMS_quality'] = df['quality']
	df['decoy']			= isDecoy
	return df[df.nXLinks == nXL]


target 		= setup(targetPath, False)
decoy 		= setup(decoyPath, True)
merged 		= concat([target, decoy], ignore_index=True)
pyProphIn	= to(targetPath, "pyProph.csv")
project 	= base(targetPath)


if dual:
	outputCols = ['transition_group_id', 'decoy', 
			'main_var_light_rt_score', 'var_heavy_rt_score', 'var_double_int_score', 'var_rel_mz_diff', 'var_openMS_quality', 
			'rt', 'mz', 'charge', 'nXLinks', 'intensity', 'width']
else:
	outputCols = ['transition_group_id', 'decoy', 
			'main_var_light_rt_score', 'var_single_int_score', 'var_rel_mz_diff', 'var_openMS_quality', 
			'rt', 'mz', 'charge', 'nXLinks', 'intensity', 'width']

merged.to_csv(
	pyProphIn, 
	index=False, 
	sep='\t',
	cols=outputCols 
	)

args = ""
for k,v in conf.iteritems():
	if v == "True":
		args += " --%s" % k
	else:
		args += " --%s=%s" % (k, v)

args += " %s " % (pyProphIn)

if dryrun:
	print "    DEBUG      "
	print "   =======     "
	print " wanted to run:"
	print "pyprophet", args
else:
	os.system("/home/johant/bin/pyprophet/bin/pyprophet " + args)
