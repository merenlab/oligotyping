# -*- coding: utf-8

# Copyright (C) 2012 - 2013, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

import os
import time

import Oligotyping.lib.fastalib as u
import Oligotyping.lib.b6lib as b6lib

from Oligotyping.utils.utils import append_file 
from Oligotyping.utils.utils import run_command
from Oligotyping.utils.utils import Multiprocessing
from Oligotyping.utils.utils import split_fasta_file
from Oligotyping.utils.utils import is_program_exist
from Oligotyping.utils.utils import check_command_output
from Oligotyping.utils.utils import get_temporary_file_name


version_error_text = '''\n
            Certain steps of decomposition requires fast searching, and it seems NCBI's BLAST tools
            installed on your system do not have the expected version number.
            
            Please make sure you have BLAST tools version 2.2.26 or higher. You can download BLAST
            tools from this address:
            
                ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
            
            If this error message does not make any sense to you, please contact your system
            administrator.
            
            Exiting.\n'''

missing_binary_error_text = '''\n
            Certain steps of decomposition requires fast searching, and it seems NCBI's BLAST tools
            are not available on your system.
            
            Please make sure 'blastn' and 'makeblastdb' binaries are in your path. You can download
            BLAST tools from here:
            
                ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
            
            If this error message does not make any sense to you, please contact your system
            administrator.
            
            Exiting.\n'''

class ModuleBinaryError(Exception):
    def __init__(self, e = None):
        Exception.__init__(self)
        self.e = e
        return
    def __str__(self):
        return 'Module Binary Error: %s' % self.e


class ModuleVersionError(Exception):
    def __init__(self, e = None):
        Exception.__init__(self)
        self.e = e
        return
    def __str__(self):
        return 'Module Version Error: %s' % self.e


class LocalBLAST:
    def __init__(self, input_fasta, target, output = None, binary = "blastn", makeblastdb = "makeblastdb", log = "/dev/null"):
        self.binary = binary
        self.params = ''
        self.input = input_fasta
        self.target = target
        self.output = output
        self.makeblastdb = makeblastdb
        self.log = log
        self.outfmt = "'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'"
        self.debug = False

        self.binary_check()        
        self.version_check()

        if not output:
            self.output = get_temporary_file_name()
        else:
            self.output = output
        
        if self.debug:
            self.search_cmd_tmpl = "%(binary)s -query %(input)s -db %(target)s -out %(output)s -outfmt %(outfmt)s %(params)s"
            self.makeblastdb_cmd_tmpl = "%(makeblastdb)s -in %(target)s -dbtype nucl"
        else:
            self.search_cmd_tmpl = "%(binary)s -query %(input)s -db %(target)s -out %(output)s -outfmt %(outfmt)s %(params)s &> %(log)s"
            self.makeblastdb_cmd_tmpl = "%(makeblastdb)s -in %(target)s -dbtype nucl &> %(log)s"
        
        self.results_dict = {}


    def get_cmd_line_params_dict(self):
        cmd_line_params_dict = {'binary': self.binary,
                                'input' : self.input,
                                'output': self.output,
                                'target': self.target,
                                'params': self.params,
                                'outfmt': self.outfmt,
                                'makeblastdb': self.makeblastdb,
                                'log': self.log}

        return cmd_line_params_dict


    def binary_check(self):
        if (not is_program_exist(self.binary)) or (not is_program_exist(self.makeblastdb)):
            raise ModuleBinaryError, missing_binary_error_text


    def version_check(self):
        version_text = check_command_output('%(binary)s -version' % self.get_cmd_line_params_dict())
        # we expect to see an output like this:
        #
        #    blastn: 2.2.26+
        #    Package: blast 2.2.26, build Feb 10 2012 02:04:27
        #
        major_blastn_version = version_text.strip().split()[1].split('.')[0]
        
        if major_blastn_version != '2':
            raise ModuleVersionError, version_error_text


    def search_parallel(self, num_processes, num_reads_per_process = 2000):
        def worker(search_cmd):
            run_command(search_cmd)
        
        mp = Multiprocessing(worker)

        input_file_parts = split_fasta_file(self.input,
                                            os.path.dirname(self.input),
                                            num_reads_per_file = num_reads_per_process)
        processes_to_run = []
        output_file_parts = []
        
        for input_file_part in input_file_parts:
            cmd_line_params_dict = self.get_cmd_line_params_dict()
            cmd_line_params_dict['input'] = input_file_part
            output_file_part = input_file_part + '.b6'
            cmd_line_params_dict['output'] = output_file_part
            output_file_parts.append(output_file_part)
            processes_to_run.append(self.search_cmd_tmpl % cmd_line_params_dict)

        self.search_cmd = processes_to_run[0]

        while 1:
            running_processes = len([p for p in mp.processes if p.is_alive()])
            
            if running_processes < num_processes and processes_to_run:
                mp.run((processes_to_run.pop(),))

            if not running_processes and not processes_to_run:
                # let the blastn program finish writing all output files.
                # FIXME: this is ridiculous. find a better solution.
                time.sleep(5)
                break

            time.sleep(1)
        
        if os.path.exists(self.output):
            os.remove(self.output)
        
        for output_file_part in output_file_parts:
            append_file(self.output, output_file_part)

        for input_file_part in input_file_parts:
            os.remove(input_file_part)
        

    def search(self, num_processes = None):
        self.search_cmd = self.search_cmd_tmpl % self.get_cmd_line_params_dict()
        run_command(self.search_cmd)


    def make_blast_db(self):
        self.makeblastdb_cmd = self.makeblastdb_cmd_tmpl % self.get_cmd_line_params_dict()
        run_command(self.makeblastdb_cmd)


    def get_results_dict(self, mismatches = None, gaps = None, min_identity = None, max_identity = None):
        results_dict = {}
        
        b6 = b6lib.B6Source(self.output)
        
        ids_with_hits = set()
        while b6.next():
            if b6.query_id == b6.subject_id:
                continue

            # following correction is to take secret gaps into consideration.
            # because we are working with reads that are supposed to be almost the
            # same length, we want query and target to be aligned 100%. sometimes it
            # is not the case, and mismatches are being calculated by the aligned
            # part of query or target. for instance if query is this:
            #
            #    ATCGATCG
            #
            # and target is this:
            #
            #   TATCGATCG
            #
            # the alignment discards the T at the beginning and gives 0 mismatches.
            # here we introduce those gaps back:
            additional_gaps = 0
            if b6.q_start != 1 or b6.s_start != 1:
                additional_gaps += (b6.q_start - 1) if b6.q_start > b6.s_start else (b6.s_start - 1)
            if additional_gaps != b6.q_len or b6.s_end != b6.s_len:
                additional_gaps += (b6.q_len - b6.q_end) if (b6.q_len - b6.q_end) > (b6.s_len - b6.s_end) else (b6.s_len - b6.s_end)

            identity_penalty = additional_gaps * 100.0 / (b6.q_len + additional_gaps)
            
            if identity_penalty:
                b6.gaps += additional_gaps
                b6.identity -= identity_penalty
                
            # done correcting the hit. carry on.

            if max_identity is not None:
                if round(b6.identity, 1) >= round(max_identity, 1):
                    continue

            if min_identity is not None:
                if round(b6.identity, 1) < round(min_identity, 1):
                    continue
            
            if mismatches is not None:
                if b6.mismatches != mismatches:
                    continue
                
            if gaps is not None:
                if b6.gaps != gaps:
                    continue
            
            # if it made here, we are interested in this one, after all.
            if b6.query_id not in ids_with_hits:
                results_dict[b6.query_id] = set()
                ids_with_hits.add(b6.query_id)
                
            results_dict[b6.query_id].add(b6.subject_id)
                
        b6.close()
        
        return results_dict
            
if __name__ == "__main__":
    try:
        u = LocalBLAST(None, None, None, None)
    except ModuleVersionError:
        raise ModuleVersionError, version_error_text
    