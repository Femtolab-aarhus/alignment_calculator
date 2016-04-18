#!/usr/bin/python3

import clustercomputing
import argparse
import time

parser = argparse.ArgumentParser(description='Start a job server');
parser.add_argument('port',type=int,help='Listen on this port');
parser.add_argument('authkey',type=str,help='Password for server');
parser.add_argument('--timeout',type=int,nargs=1,default=[120],help='TODO');

args = parser.parse_args();

address = ('', args.port);
authkey = bytes(args.authkey,'utf-8');
timeout = args.timeout[0];

manager,srv = clustercomputing.start_server(address,authkey,timeout);

print("Server started. Press any key to exit.");
input();
manager.shutdown();




