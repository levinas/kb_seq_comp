#BEGIN_HEADER
import os
import sys
import traceback
import json
import logging
import subprocess
import tempfile
import uuid
import hashlib

from datetime import datetime

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biokbase.workspace.client import Workspace as workspaceService


logging.basicConfig(format="[%(asctime)s %(levelname)s %(name)s] %(message)s",
                    level=logging.DEBUG)
logger = logging.getLogger(__name__)

#END_HEADER


class SeqComparison:
    '''
    Module Name:
    SeqComparison

    Module Description:
    A KBase module: SeqComparison
This sample module contains one small method - filter_contigs.
    '''

    ######## WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    #BEGIN_CLASS_HEADER
    workspaceURL = None

    # target is a list for collecting log messages
    def log(self, target, message):
        # we should do something better here...
        if target is not None:
            target.append(message)
        logger.info(message)

    def contigset_to_fasta(self, contigset, fasta_file):
        records = []
        for contig in contigset['contigs']:
            record = SeqRecord(Seq(contig['sequence']), id=contig['id'], description='')
            records.append(record)
        SeqIO.write(records, fasta_file, "fasta")

    def create_temp_json(self, attrs):
        f = tempfile.NamedTemporaryFile(delete=False)
        outjson = f.name
        f.write(json.dumps(attrs))
        f.close()
        return outjson
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.workspaceURL = config['workspace-url']
        self.scratch = os.path.abspath(config['scratch'])
        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)
        #END_CONSTRUCTOR
        pass

    def run_dnadiff(self, ctx, params):
        # ctx is the context object
        # return variables are: output
        #BEGIN run_dnadiff

        logger.info("params={}".format(json.dumps(params)))

        token = ctx["token"]
        ws = workspaceService(self.workspaceURL, token=token)
        wsid = None

        genomeset = None
        if "input_genomeset_ref" in params and params["input_genomeset_ref"] is not None:
            logger.info("Loading GenomeSet object from workspace")
            objects = ws.get_objects([{"ref": params["input_genomeset_ref"]}])
            genomeset = objects[0]["data"]
            wsid = objects[0]['info'][6]

        genome_refs = []
        if genomeset is not None:
            for param_key in genomeset["elements"]:
                genome_refs.append(genomeset["elements"][param_key]["ref"])
            logger.info("Genome references from genome set: {}".format(genome_refs))

        if "input_genome_refs" in params and params["input_genome_refs"] is not None:
            for genome_ref in params["input_genome_refs"]:
                if genome_ref is not None:
                    genome_refs.append(genome_ref)

        logger.info("Final list of genome references: {}".format(genome_refs))
        if len(genome_refs) < 2:
            raise ValueError("Number of genomes should be more than 1")
        if len(genome_refs) > 10:
            raise ValueError("Number of genomes exceeds 10, which is too many for dnadiff")

        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()*1000)
        output_dir = os.path.join(self.scratch, 'output.'+str(timestamp))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        genome_names = []
        fasta_files = []
        for pos, ref in enumerate(genome_refs):
            logger.info("Loading Genome object from workspace for ref: ".format(ref))

            obj = ws.get_objects([{"ref": ref}])[0]
            data = obj["data"]
            info = obj["info"]
            wsid = wsid or info[6]
            type_name = info[2].split('.')[1].split('-')[0]
            logger.info("type_name = {}".format(type_name))

            # if KBaseGenomes.ContigSet
            if type_name == 'Genome':
                # logger.debug("genome = {}".format(json.dumps(data)))
                genome_names.append(data.get("scientific_name", "") + " ({})".format(ref))
                contigset_ref = data["contigset_ref"]
                obj = ws.get_objects([{"ref": contigset_ref}])[0]
                data = obj["data"]
                info = obj["info"]
                # logger.debug("data = {}".format(json.dumps(data)))
            else:
                genome_names.append(ref.split('/')[1] + " ({})".format(ref.split('/')[0]))


            fasta_name = os.path.join(output_dir, "{}.fa".format(pos+1))
            self.contigset_to_fasta(data, fasta_name)
            fasta_files.append(fasta_name)

            # data_ref = str(info[6]) + "/" + str(info[0]) + "/" + str(info[4])

            # logger.info("info = {}".format(json.dumps(info)))
            # logger.info("data = <<<<<<{}>>>>>>".format(json.dumps(data)))

        logger.info("fasta_files = {}".format(fasta_files))

        logger.info("Run DNAdiff:")

        out_json = os.path.join(output_dir, 'out.json')
        cmd = ['dnadiff_genomes', '-o', 'out', '-j', out_json]
        cmd += fasta_files

        # logger.info("CMD: {}".format(' '.join(cmd)))
        p = subprocess.Popen(cmd,
                             cwd = self.scratch,
                             stdout = subprocess.PIPE,
                             stderr = subprocess.STDOUT, shell = False)

        console = []
        while True:
            line = p.stdout.readline()
            if not line: break
            self.log(console, line.replace('\n', ''))

        p.stdout.close()
        p.wait()
        logger.debug('return code: {}'.format(p.returncode))
        if p.returncode != 0:
            raise ValueError('Error running dnadiff, return code: {}\n\n{}'.format(p.returncode, '\n'.join(console)))

        logger.info('Genomes/ContigSets compared with DNAdiff:\n')
        for pos, name in enumerate(genome_names):
            logger.info('  {}: {}\n'.format(pos+1, name))

        with open(out_json) as comp_file:
            comparisons = json.load(comp_file)

        comp_data = { 'genome_names': genome_names, 'genome_comparisons': comparisons }


        # provenance
        input_ws_objects = []
        if "input_genomeset_ref" in params and params["input_genomeset_ref"] is not None:
            input_ws_objects.append(params["input_genomeset_ref"])
        if "input_genome_refs" in params and params["input_genome_refs"] is not None:
            for genome_ref in params["input_genome_refs"]:
                if genome_ref is not None:
                    input_ws_objects.append(genome_ref)

        provenance = None
        if "provenance" in ctx:
            provenance = ctx["provenance"]
        else:
            logger.info("Creating provenance data")
            provenance = [{"service": "SeqComparison",
                           "method": "run_dnadiff",
                           "method_params": [params]}]

        provenance[0]["input_ws_objects"] = input_ws_objects
        provenance[0]["description"] = "Sequence comparison using DNAdiff"


        # save the report object
        comp_obj_info = ws.save_objects({
            'id': wsid, # set the output workspace ID
            'objects':[{'type': 'ComparativeGenomics.SeqCompOutput',
                        'data': comp_data,
                        'name': params['output_report_name'],
                        'meta': {},
                        'provenance': provenance}]})


        # reportObj = {
        #     'objects_created':[{'ref':params['workspace_name']+'/'+params['output_report_name'], 'description':'DNAdiff report'}],
        #     'text_message': report
        # }

        # reportName = '{}.report.{}'.format('run_dnadiff', hex(uuid.getnode()))
        # report_obj_info = ws.save_objects({
        #         # 'workspace': params["workspace_name"],
        #     'id': wsid,
        #     'objects': [
        #         {
        #             'type': 'KBaseReport.Report',
        #             'data': reportObj,
        #             'name': reportName,
        #             'meta': {},
        #             'hidden': 1,
        #             'provenance': provenance
        #         }
        #     ]})[0]


        # shutil.rmtree(output_dir)

        # output = {"report_name": reportName, 'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]) }
        output = {}

        #END run_dnadiff

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_dnadiff return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
