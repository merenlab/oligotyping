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
import copy
import io

import Oligotyping.lib.fastalib as u
import Oligotyping.lib.b6lib as b6lib

from Oligotyping.utils.utils import append_file 
from Oligotyping.utils.utils import run_command
from Oligotyping.utils.utils import Multiprocessing
from Oligotyping.utils.utils import split_fasta_file
from Oligotyping.utils.utils import is_program_exist
from Oligotyping.utils.utils import check_command_output
from Oligotyping.utils.utils import get_temporary_file_name
from Oligotyping.utils.utils import remove_white_space_mask_from_B6_entry
from Oligotyping.utils.aligner import nw_align 


biopython_error_text = '''\n
            You need 'BioPython' module in your Python path to run this software.

            You can get more information about the installation here:

                http://biopython.org/wiki/Biopython

            Exiting.\n'''

version_error_text = '''\n
            Certain steps of decomposition require fast searching, and it seems NCBI's BLAST tools
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


class MissingModuleError(Exception):
    def __init__(self, e = None):
        Exception.__init__(self)
        self.e = e
        return
    def __str__(self):
        return 'Missing Module Error: %s' % self.e


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


try:
    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML
except:
    raise MissingModuleError(biopython_error_text)


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

        self.binary_check()        
        self.version_check()

        if not output:
            self.output = get_temporary_file_name()
        else:
            self.output = output
        
        self.search_cmd_tmpl = "%(binary)s -query %(input)s -db %(target)s -out %(output)s -outfmt %(outfmt)s %(params)s >> %(log)s 2>&1"
        self.makeblastdb_cmd_tmpl = "%(makeblastdb)s -in %(target)s -dbtype nucl >> %(log)s 2>&1"
        
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
            raise ModuleBinaryError(missing_binary_error_text)


    def version_check(self):
        version_text = check_command_output('%(binary)s -version' % self.get_cmd_line_params_dict()).decode("utf-8") 
        # we expect to see an output like this:
        #
        #    blastn: 2.2.26+
        #    Package: blast 2.2.26, build Feb 10 2012 02:04:27
        #
        major_blastn_version = version_text.strip().split()[1].split('.')[0]
        
        if major_blastn_version != '2':
            raise ModuleVersionError(version_error_text)


    def search_parallel(self, num_processes, num_reads_per_process = 2000, keep_parts = False):
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

            time.sleep(5)
        
        if os.path.exists(self.output):
            os.remove(self.output)
        
        for output_file_part in output_file_parts:
            if keep_parts:
                append_file(self.output, output_file_part, remove_source = False)
            else:
                append_file(self.output, output_file_part)

        if not keep_parts:
            for input_file_part in input_file_parts:
                os.remove(input_file_part)
        

    def search(self, num_processes = None):
        self.search_cmd = self.search_cmd_tmpl % self.get_cmd_line_params_dict()
        run_command(self.search_cmd)


    def make_blast_db(self):
        self.makeblastdb_cmd = self.makeblastdb_cmd_tmpl % self.get_cmd_line_params_dict()
        run_command(self.makeblastdb_cmd)


    def get_results_dict(self, mismatches = None, gaps = None, min_identity = None, max_identity = None, penalty_for_terminal_gaps = True):
        results_dict = {}
        
        b6 = b6lib.B6Source(self.output)

        ids_with_hits = set()
        while b6.next():
            if b6.entry.query_id == b6.entry.subject_id:
                continue

            if penalty_for_terminal_gaps:
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
                if b6.entry.q_start != 1 or b6.entry.s_start != 1:
                    additional_gaps += (b6.entry.q_start - 1) if b6.entry.q_start > b6.entry.s_start else (b6.entry.s_start - 1)
                if additional_gaps != b6.entry.q_len or b6.entry.s_end != b6.entry.s_len:
                    additional_gaps += (b6.entry.q_len - b6.entry.q_end) if (b6.entry.q_len - b6.entry.q_end) > (b6.entry.s_len - b6.entry.s_end) else (b6.entry.s_len - b6.entry.s_end)
    
                identity_penalty = additional_gaps * 100.0 / (b6.entry.q_len + additional_gaps)
                
                if identity_penalty:
                    b6.entry.gaps += additional_gaps
                    b6.entry.identity -= identity_penalty
                    
                # done correcting the hit. carry on.

            if max_identity is not None:
                if round(b6.entry.identity, 1) >= round(max_identity, 1):
                    continue

            if min_identity is not None:
                if round(b6.entry.identity, 1) < round(min_identity, 1):
                    continue
            
            if mismatches is not None:
                if b6.entry.mismatches != mismatches:
                    continue
                
            if gaps is not None:
                if b6.entry.gaps != gaps:
                    continue
            
            # if it made here, we are interested in this one, after all.
            if b6.entry.query_id not in ids_with_hits:
                results_dict[b6.entry.query_id] = set()
                ids_with_hits.add(b6.entry.query_id)
                
            results_dict[b6.entry.query_id].add(b6.entry.subject_id)
                
        b6.close()
        
        return results_dict


    def get_fancy_results_dict(self, max_per_query = 10, defline_white_space_mask = None):
        b6 = b6lib.B6Source(self.output)

        input_fasta = u.SequenceSource(self.input)
        target_db = u.SequenceSource(self.target)

        query_counts = {}
        fancy_results_dict = {}

        while next(b6):
            if b6.entry.query_id not in query_counts:
                query_counts[b6.entry.query_id] = 1

            if query_counts[b6.entry.query_id] - 1 == max_per_query:
                continue
            else:
                query_counts[b6.entry.query_id] += 1

            if b6.entry.query_id not in fancy_results_dict:
                fancy_results_dict[b6.entry.query_id] = []

            query_seq = input_fasta.get_seq_by_read_id(b6.entry.query_id).replace('-', '')
            target_seq = target_db.get_seq_by_read_id(b6.entry.subject_id)
            
            if defline_white_space_mask:
                b6.entry = remove_white_space_mask_from_B6_entry(b6.entry, defline_white_space_mask)
            
            # parts that were aligned during the search are being aligned to each other to generate
            # hsp_match data to include into results
            query_aligned, target_aligned = nw_align(query_seq[int(b6.entry.q_start) - 1:int(b6.entry.q_end)],\
                                                         target_seq[int(b6.entry.s_start) - 1:int(b6.entry.s_end)]) 

            query_aligned, target_aligned = query_aligned.upper(), target_aligned.upper()

            coverage = (b6.entry.q_end - (b6.entry.q_start - 1)) * 100.0 / b6.entry.q_len
            hsp_match = ''.join(['|' if query_aligned[i] == target_aligned[i] else ' ' for i in range(0, len(query_aligned))])
            
            entry = copy.deepcopy(b6.entry)
            entry.coverage = coverage
            entry.hsp_query = query_aligned
            entry.hsp_subject = target_aligned
            entry.hsp_match = hsp_match
            
            entry = remove_white_space_mask_from_B6_entry(entry)

            fancy_results_dict[entry.query_id].append(entry)

        return fancy_results_dict


class RemoteBLAST:
    def __init__(self):
        pass

    def search(self, sequence, output_file = None):
        result_handle = NCBIWWW.qblast("blastn", "nt", sequence.replace('-', ''))
        result = result_handle.read()

        if output_file:
            open(output_file, "w").write(result)
 
        return io.StringIO(result)


    def get_fancy_results_list(self, blast_results, num_results = 20):
        blast_results_list = []
    
        blast_record = list(NCBIXML.parse(blast_results))[0]
        num_results = len(blast_record.alignments) if len(blast_record.alignments) < num_results else num_results
    
        for i in range(0, num_results):
            entry = b6lib.B6Entry()
            entry.q_len = int(blast_record.query_length)
            entry.query_length = entry.q_len
            
            alignment = blast_record.alignments[i]
            hsp = alignment.hsps[0]
         
            entry.hit_def = alignment.hit_def   
            entry.subject_id = entry.hit_def
            entry.accession = alignment.accession
            entry.ncbi_link = 'http://www.ncbi.nlm.nih.gov/nuccore/%s' % entry.accession
            entry.hsp_query = hsp.query
            entry.hsp_match = hsp.match
            entry.hsp_subject = hsp.sbjct
    
            entry.identity = len([x for x in hsp.match if x == '|']) * 100.0 / len(entry.hsp_query)
            entry.coverage = len(hsp.query) * 100.0 / entry.query_length
    
            blast_results_list.append(entry)
    
        try:
            blast_results.close()
        except:
            pass
    
        return blast_results_list


if __name__ == "__main__":
    try:
        u = LocalBLAST(None, None)
    except ModuleVersionError:
        raise ModuleVersionError(version_error_text)
    
