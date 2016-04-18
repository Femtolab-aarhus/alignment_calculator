import multiprocessing
import multiprocessing.managers
import queue
import time
import socket

class NoMoreJobs(Exception):
    pass;

class jobserver(object):
    def __init__(s,timeout=60):
        s.job_queue = multiprocessing.Queue();
        s.result_queue = multiprocessing.Queue();
        s.curr_jobs = dict();
        s.lock = multiprocessing.Lock();
        s.timeout = timeout;
        s.id_counter = 0;
        s.shutdown = False;
        
    def close(s):
        s.shutdown = True;
        
    def _vacuum(s):
        now = time.time();
        timedout = [];
        for job_id,package in s.curr_jobs.items():
            timestamp, job = package;
            if (now - timestamp > s.timeout):
                timedout.append((job_id,job));
        for job_id,job in timedout:
            del s.curr_jobs[job_id];
            s.job_queue.put((job_id,job));

    def submit_job(s,job):
        if (s.shutdown):
            return None;
        with s.lock:
            job_id = s.id_counter;
            s.id_counter += 1;
        s.job_queue.put((job_id,job));
        return job_id;


    def get_job(s,timeout):
        try:
            job_id,job = s.job_queue.get(timeout=timeout);
        except queue.Empty:
            with s.lock:
                s._vacuum();
                try:
                    job_id,job = s.job_queue.get_nowait();
                except queue.Empty:
                    if (s.shutdown and len(s.curr_jobs) == 0):
                        raise NoMoreJobs();
                    else:
                        return None;
        with s.lock:
            s.curr_jobs[job_id] = (time.time(),job);
        return job_id,job;
        
    def submit_result(s,job_id,result,credits=None):
        with s.lock:
            try:
                del s.curr_jobs[job_id];
            except KeyError:
                print("Unknown job result submitted");
                return;
            if (credits):
                print(credits);
        s.result_queue.put((job_id,result));
        
    def get_result(s,timeout=None):
        return s.result_queue.get(timeout=timeout);


class jobserver_manager(multiprocessing.managers.BaseManager):
            pass


def start_server(address,authkey,timeout=60):
    srv = jobserver(timeout);
    jobserver_manager.register('jobserver',callable=lambda:srv);
    manager = jobserver_manager(address=address,authkey=authkey);
    manager.start();
    return manager,srv;

def connect_to_sever(address,authkey):
    jobserver_manager.register('jobserver');
    manager = jobserver_manager(address=address,authkey=authkey);
    manager.connect();
    srv = manager.jobserver();
    return manager,srv;
    

