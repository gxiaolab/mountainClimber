#!/user/bin/python -tt

import os
import sys
import re
from subprocess import Popen, PIPE


def run_command(cmd, stdin=0, stdoutfile=0, printflag=1):
	"""run a subprocess"""
	if printflag == 1:
		print ' '.join(cmd)
	if stdin == 0 and stdoutfile == 0:
		p = Popen(cmd, stdout=PIPE, stderr=PIPE)
		stdout, stderr = p.communicate()
	elif stdoutfile == 0:
		p = Popen(cmd, stdout=PIPE, stderr=PIPE, stdin=PIPE)
		stdout, stderr = p.communicate(stdin)
	else:
		p = Popen(cmd, stdout=stdoutfile, stderr=PIPE)
		pcode = p.wait()
		stdoutfile.flush()
		stdout = 'NA'
		stderr = 'NA'

	if p.returncode != 0:
		if stdoutfile == 0:
			sys.stderr.write('command failed')
			sys.stderr.write(stdout)
			sys.stderr.write(stderr)
			sys.exit(1)
		else:
			sys.stderr.write('command failed')
			sys.stderr.write(stderr)
			sys.exit(1)

	return stdout, stderr


def sort_bedfile(infile, outfile):
	"""sort bed file"""
	o = open(outfile, 'w')
	o.write('track name=\"' + os.path.basename(infile).replace('.bed', '') + '\"\n')
	cmd = ['sort', '-k1,1', '-k2,2n', infile]
	stdout, stderr = run_command(cmd)
	o.write(stdout)
	o.close()
