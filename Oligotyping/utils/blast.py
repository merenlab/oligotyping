# -*- coding: utf-8

# Copyright (C) 2012 - 2013, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

import Oligotyping.lib.b6lib as b6lib

from Oligotyping.utils.utils import run_command
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
    def __init__(self, input_fasta, target, params, output = None, binary = "blastn", makeblastdb = "makeblastdb", log = "/dev/null"):
        self.binary = binary
        self.params = params
        self.input = input_fasta
        self.target = target
        self.output = output
        self.makeblastdb = makeblastdb
        self.log = log

        self.cmd_line_params_dict = {'binary': self.binary,
                                     'input' : self.input,
                                     'output': self.output,
                                     'target': self.target,
                                     'params': self.params,
                                     'makeblastdb': self.makeblastdb,
                                     'log': self.log}

        self.binary_check()        
        self.version_check()

        if not output:
            self.output = get_temporary_file_name()
        else:
            self.output = output
        
        self.search_cmd_tmpl = "%(binary)s -query %(input)s -db %(target)s -out %(output)s -outfmt 6 %(params)s &> %(log)s"
        self.makeblastdb_cmd_tmpl = "%(makeblastdb)s -in %(target)s -dbtype nucl &> %(log)s"
        
        self.results_dict = {}


    def binary_check(self):
        if (not is_program_exist(self.binary)) or (not is_program_exist(self.makeblastdb)):
            raise ModuleBinaryError, missing_binary_error_text


    def version_check(self):
        version_text = check_command_output('%(binary)s -version' % self.cmd_line_params_dict)
        # we expect to see an output like this:
        #
        #    blastn: 2.2.26+
        #    Package: blast 2.2.26, build Feb 10 2012 02:04:27
        #
        major_blastn_version = version_text.strip().split()[1].split('.')[0]
        
        if major_blastn_version != '2':
            raise ModuleVersionError, version_error_text


    def search(self):
        self.search_cmd = self.search_cmd_tmpl % self.cmd_line_params_dict
        run_command(self.search_cmd)


    def make_blast_db(self):
        self.makeblastdb_cmd = self.makeblastdb_cmd_tmpl % self.cmd_line_params_dict
        run_command(self.makeblastdb_cmd)


    def get_results_dict(self, mismatches = 1, gaps = 0):
        results_dict = {}
        
        b6 = b6lib.B6Source(self.output)
        
        ids_with_hits = set()
        while b6.next():
            if b6.query_id == b6.subject_id:
                continue
            if b6.mismatches == mismatches and b6.gaps == gaps:
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
    