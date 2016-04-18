#!/usr/bin/python3

import propagation
import interaction
import numpy 
from numpy import pi,exp,sqrt,log,sin,cos
import argparse
import time
import datetime
import multiprocessing
import clustercomputing
import socket

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Work client for prepare_thermal');

    parser.add_argument('host',type=str,help='Server ip or host name');
    parser.add_argument('authkey',type=str,help='Password for job server');
    parser.add_argument('--port',nargs=1,default=[2000],type=int,help='Server port number');
        
    args = parser.parse_args();

    host = args.host;
    port = args.port[0];

    address = (host,port);
    authkey = bytes(args.authkey,'utf-8');

    manager,job_server = clustercomputing.connect_to_sever(address,authkey);

    pool = multiprocessing.Pool();

    my_name = socket.getfqdn();

    jobs = [];

    job_id, args = job_server.get_job(1);
    args.append(pool);
    async_results = propagation.compute_transfer_matrix_async(*args);
    jobs.append((job_id,async_results));

    try:
        while True:
            before = time.time();
            try:
                job = job_server.get_job(10);
                if (job):
                    job_id, args = job;
                else:
                    print("Job server didn't serve any jobs.");
                    break;
            except clustercomputing.NoMoreJobs:
                print("Job server closed.");
                break;

            args.append(pool);
            async_results = propagation.compute_transfer_matrix_async(*args);
            jobs.append((job_id,async_results));

            job_id,async_results = jobs.pop(0);
        
            transfer_matrix = propagation.collect_async_transfer_matrix(async_results);
            job_server.submit_result(job_id,transfer_matrix,my_name);
            print(datetime.timedelta(seconds=time.time()-before));

    except KeyboardInterrupt:
        print("Finishing remaining jobs");

 
    for job_id,async_results in jobs:
        transfer_matrix = propagation.collect_async_transfer_matrix(async_results);
        job_server.submit_result(job_id,transfer_matrix,my_name);
        print("Finished remaining job.");

    pool.close();
    pool.join();


