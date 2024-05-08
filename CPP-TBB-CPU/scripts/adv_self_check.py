#!/usr/bin/env python
#
# Copyright (C) 2016 Intel Corporation
#
# This software and the related documents are Intel copyrighted materials, and your use of them
# is governed by the express license under which they were provided to you ("License"). Unless
# the License provides otherwise, you may not use, modify, copy, publish, distribute, disclose
# or transmit this software or the related documents without Intel's prior written permission.
#
# This software and the related documents are provided as is, with no express or implied
# warranties, other than those that are expressly stated in the License.
#

import sys
import json
import os
import shlex, subprocess
import csv

from tempfile import gettempdir
from optparse import OptionParser

from self_check_common import *

CL_TOOL = "advisor"
PRODUCT_ABBR = "advisor"
INTERNAL_PREFIX = "advixe-"

ICL_CPU_EXE_NAME  = "o3memog"
ICL_CPU_EXE_OPTS  = shlex.split("-i ../bin/data/BrionUH -W 10 -B 20 -O 0.8 -n 8 -l 5 -r 0.1 -s 32 -L 4")
PROJECT_NAME = 'gcl_cpu_result'

DPCPP_CPU_EXE_NAME  = "mtx.mult.dpcpp.cpu"
DPCPP_CPU_EXE_OPTS  = ""

DPCPP_GPU_EXE_NAME  = "mtx.mult.dpcpp.gpu"
DPCPP_GPU_EXE_OPTS  = ""

# Analysis types
ICL_CPU_SURVEY_ANALYSIS = "Survey analysis"
ICL_CPU_TRIPCOUNTS_ANALYSIS = "Tripcounts analysis"
ICL_CPU_MAP_ANALYSIS = "Map analysis"
ICL_CPU_DEPENDENCIES_ANALYSIS = "Dependencies analysis"
ICL_CPU_PROJECTION_ANALYSIS = "Projection analysis"
ICL_CPU_ROOFLINE_ANALYSIS = "Roofline analysis"

DPCPP_CPU_SURVEY_ANALYSIS = "Survey analysis"
DPCPP_CPU_TRIPCOUNTS_FLOP_ANALYSIS = "Tripcounts/FLOP analysis"

DPCPP_GPU_OCL_ROOFLINE_ANALYSIS = "OpenCL Roofline analysis"
DPCPP_GPU_L0_ROOFLINE_ANALYSIS = "Level-zero Roofline analysis"

COLLECTION = "Collection"
REPORT = "Report"

COPYRIGHTS_ADV = "Intel(R) Advisor Self Check Utility"

FOLDER_NOT_AVAILABLE = "Input directory for log is not available. Use an existing directory."
SYSTEM_READY = "The system is ready to be used for performance analysis with Intel(R) Advisor."
SYSTEM_NOT_READY = "The check observed a product failure on your system.\n"
SYSTEM_NOT_READY += "Review errors in the output above to fix a problem or contact Intel technical support."
UNKNOWN_FAIL = "An unknown failure has occurred. Attempt the action again or contact Intel technical support."
SEE_WARNINGS = "Review warnings in the output above to find product limitations, if any."

HELP_LOG_DIR = "path to directory where to store log"

ERROR_TAG = PRODUCT_ABBR + ": Error:"
WARNING_TAG = PRODUCT_ABBR + ": Warning:"
CANNOT_LOCATE_FILE_WARNING = "Cannot locate file"
CANNOT_LOCATE_DEBUGGING_INFORMATION_WARNING = "Cannot locate debugging information for file"


advisor_icl_cpu_test_descriptor = [
    {
        'name': ICL_CPU_SURVEY_ANALYSIS,
        'analysis': ['survey',],
        'analysis_opts': '',
        'report': ['survey',],
    },
    {
        'name': ICL_CPU_TRIPCOUNTS_ANALYSIS,
        'analysis': ['tripcounts'],
        'analysis_opts': ['--flop', '--stacks', '--no-auto-finalize'],
        'report': '',
    },
    {
        'name': ICL_CPU_MAP_ANALYSIS,
        'analysis': ['map',],
        'analysis_opts': '',
        'report': ['map',],
    },
    {
        'name': ICL_CPU_DEPENDENCIES_ANALYSIS,
        'analysis': ['dependencies',],
        'analysis_opts': ['--filter-reductions'],
        'report': ['dependencies',],
    },
    # {
    #     'name': ICL_CPU_PROJECTION_ANALYSIS,
    #     'analysis': ['projection'],
    #     'analysis_opts': ['--config=gen9_gt2'],
    #     'report': '',
    # },
    # {
    #     'name': ICL_CPU_ROOFLINE_ANALYSIS,
    #     'analysis': ['roofline'],
    #     'analysis_opts': ['--stacks'],
    #     'report': ['roofline'],
    # },
]

advisor_dpcpp_cpu_test_descriptor = [
    {
        'name': DPCPP_CPU_SURVEY_ANALYSIS,
        'analysis': ['survey',],
        'analysis_opts': '',
        'report': '',
    },
    {
        'name': DPCPP_CPU_TRIPCOUNTS_FLOP_ANALYSIS,
        'analysis': ['tripcounts',],
        'analysis_opts': ['--flop', '--stacks'],
        'report': ['survey',],
    },
]

advisor_dpcpp_gpu_test_descriptor = [
    {
        'name': DPCPP_GPU_OCL_ROOFLINE_ANALYSIS,
        'analysis': ['roofline',],
        'analysis_opts': ['--profile-gpu'],
        'report': ['roofline'],
    },
    {
        'name': DPCPP_GPU_L0_ROOFLINE_ANALYSIS,
        'analysis': ['roofline',],
        'analysis_opts': ['--profile-gpu'],
        'report': ['roofline'],
    },    
]


class AdvisorIclCPUTest(TestCase):
    def run(self):
        self.status = 'OK'
        self.log.to_stdout("=" * 80)
        self.log.to_stdout("Running ICL CPU test:")
        os.environ["OMP_NUM_THREADS"] = "8"

        for test_descriptor in advisor_icl_cpu_test_descriptor:
            self.log.to_log("=" * 80)
            test_title = test_descriptor['name']
            self.log.to_stdout(test_title)
            self.log.to_stdout(self.state.output_indent + 'Running collection...')

            title, status, important_messages = self.run_collection(test_descriptor)

            self.log.to_stdout(title + ': ' + status)

            for message in important_messages:
                self.log.to_stdout(message)

        self.log.to_stdout('')

        for test_descriptor in advisor_icl_cpu_test_descriptor:
            if self.status == 'OK' and test_descriptor['report'] != '':
                report_title = test_descriptor['name']
                self.log.to_log("-" * 80)
                self.log.to_stdout(report_title + ' report')
                title, status, important_messages = self.run_report(test_descriptor)
                self.log.to_stdout(title + ': ' + status)

                for message in important_messages:
                    self.log.to_stdout(message)

        # if self.status == 'OK':
        #     self.remove_result_dir(test_descriptor)
        
    def run_collection(self, test_descriptor):
        title = self.state.output_indent + COLLECTION

        args = [self.state.cl_path]
        args += ['--collect']
        args += test_descriptor['analysis']
        args += test_descriptor['analysis_opts']
        args += ['--project-dir']

        if test_descriptor['name'] == ICL_CPU_ROOFLINE_ANALYSIS:
            result_dir_abs_path = self.state.append_to_work_dir(f'{PROJECT_NAME}/rl000')
        else:
            result_dir_abs_path = self.state.append_to_work_dir(PROJECT_NAME)

        args += [result_dir_abs_path]

        # app_dir = os.path.join(self.state.product_dir, 'samples', 'en', 'C++', 'adv_check')
        app_name = ICL_CPU_EXE_NAME
        if sys.platform == 'win32':
            app_name += '.exe'
        # the app is in the folder bin/, previous to the file folder
        app_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '../bin', app_name)
        
        args += ['--', app_path]
        args += ICL_CPU_EXE_OPTS

        ret_code, stdout, stderr = self.subprocess_wrapper(args, indent=self.state.stdout_indent)

        important_messages = self.get_errors_from_stderr(stderr)

        if ret_code != 0:
            self.status = 'FAIL'
            status = 'Fail'
        else:
            status = 'Ok'

        return title, status, important_messages

    def run_report(self, test_descriptor):
        title = self.state.output_indent + REPORT

        if test_descriptor['name'] == ICL_CPU_ROOFLINE_ANALYSIS:
            result_dir_abs_path = self.state.append_to_work_dir(f'{PROJECT_NAME}/rl000')
        else:
            result_dir_abs_path = self.state.append_to_work_dir(PROJECT_NAME)

        args = [self.state.cl_path,]
        args += ['--report']

        if test_descriptor['name'] == ICL_CPU_ROOFLINE_ANALYSIS:
            args += ['survey']
        else:
            args += test_descriptor['report']

        args += ['--show-all-columns']
        args += ['--format', 'csv',]
        args += ['--project-dir']
        args += [result_dir_abs_path]

        ret_code, report, stderr = self.subprocess_wrapper(args, indent=self.state.stdout_indent)

        important_messages = self.get_errors_from_stderr(stderr)

        if ret_code != 0:
            self.status = 'FAIL'
            status = 'Fail'
            return title, status, important_messages

        # if test_descriptor['name'] == ICL_CPU_ROOFLINE_ANALYSIS:        
        #     tmp_str = ''.join(map(str, 'survey'))
        # else:
        #     tmp_str = ''.join(map(str, test_descriptor['report']))

        # if test_descriptor['report'] == ['survey']:
        #         res_folder = 'hs000'
        # elif test_descriptor['report'] == ['map']:
        #         res_folder = 'mp000'
        # elif test_descriptor['report'] == ['dependencies']:
        #         res_folder = 'dp000'
        # elif test_descriptor['report'] == ['tripcounts']:
        #         res_folder = 'trc000'
        # elif test_descriptor['report'] == ['projection']:
        #         res_folder = 'pp000'
        # elif test_descriptor['report'] == ['roofline']:
        #         res_folder = 'hs000'

        # csvReport = os.path.join(result_dir_abs_path, 'e000', res_folder, '{}-{}.csv'.format(CL_TOOL, tmp_str))

        # if test_descriptor['report'] == ['survey']:
        #     LOOP_SOURCE_LOCATION = 'mmult_omp.cpp:72'
        #     EXP_TRIP_COUNT = '512'
        #     EXP_CALL_COUNT = '262144'
        #     COMPUTE_MEMORY_METRICS = [ 'Self GINTOP',
        #                                'Self GFLOP',
        #                                'Self Memory GB', ]
        #     loopFound = False
        #     with open(csvReport, 'r') as csvFile:
        #         # Skipping first 6 irrelevant lines...
        #         for i in range(6):
        #             next(csvFile)
        #         reader = csv.DictReader(csvFile)
        #         for row in reader:
        #             if row['Source Location'] == LOOP_SOURCE_LOCATION:
        #                 loopFound = True
        #                 if row['Average Trip Count'] != EXP_TRIP_COUNT:
        #                     self.log.to_log(ICL_CPU_SURVEY_ANALYSIS + " FAILED: Incorrect \'Average Trip Count\' for C++ loop")
        #                     self.status = 'FAIL'
        #                     status = 'Fail'
        #                     return title, status, important_messages
        #                 else:
        #                     self.log.to_log(ICL_CPU_SURVEY_ANALYSIS + ": Found correct \'Average Trip Count\' for C++ loop: {}" \
        #                         .format(row['Average Trip Count']))
        #                 if row['Call Count'] != EXP_CALL_COUNT:
        #                     self.log.to_log(ICL_CPU_SURVEY_ANALYSIS + " FAILED: Incorrect \'Call Count\' for C++ loop")
        #                     self.status = 'FAIL'
        #                     status = 'Fail'
        #                     return title, status, important_messages
        #                 else:
        #                     self.log.to_log(ICL_CPU_SURVEY_ANALYSIS + ": Found correct \'Call Count\' for C++ loop: {}" \
        #                         .format(row['Call Count']))
        #                 for name in COMPUTE_MEMORY_METRICS:
        #                     if not name in row or not row[name]:
        #                         self.log.to_log((ICL_CPU_SURVEY_ANALYSIS + " FAILED: No or empty '{}'" \
        #                             + " for C++ loop").format(name))
        #                         self.status = 'FAIL'
        #                         status = 'Fail'
        #                         return title, status, important_messages
        #                     else:
        #                         self.log.to_log((ICL_CPU_SURVEY_ANALYSIS + ": Found non-empty '{}'" \
        #                             + " value for C++ loop: {}").format(name, row[name]))
        #                 break
        #     if not loopFound:
        #         self.log.to_log(ICL_CPU_SURVEY_ANALYSIS + " FAIILED: Can't find C++ loop in the Survey report")
        #         self.status = 'FAIL'
        #         status = 'Fail'
        #         return title, status, important_messages

        # elif test_descriptor['report'] == ['map']:
        #     with open(csvReport, 'r') as csvFile:
        #         MAP_Report = csvFile.read()
        #     if 'Mixed Strides' not in MAP_Report:
        #         self.log.to_log(ICL_CPU_MAP_ANALYSIS + " FAILED: Can't find expected access pattern in the MAP report")
        #         self.status = 'FAIL'
        #         status = 'Fail'
        #         return title, status, important_messages
        #     else:
        #         self.log.to_log(ICL_CPU_MAP_ANALYSIS + ": Found expected access pattern in the MAP report: 'Mixed Strides'")

        # elif test_descriptor['report'] == ['dependencies']:
        #     with open(csvReport, 'r') as csvFile:
        #         dependenciesReport = csvFile.read()
        #     if 'mmult_omp.cpp:73' not in dependenciesReport:
        #         self.log.to_log(ICL_CPU_DEPENDENCIES_ANALYSIS + " FAILED: Can't find loop site location in the Dependencies report")
        #         self.status = 'FAIL'
        #         status = 'Fail'
        #         return title, status, important_messages
        #     else:
        #         self.log.to_log(ICL_CPU_DEPENDENCIES_ANALYSIS + ": Found loop site location in the Dependencies report: "
        #             + "'mmult_omp.cpp:73' as expected ")

        # elif test_descriptor['report'] == ['roofline']:
        #     LOOP_SOURCE_LOCATION = 'mmult_omp.cpp:72'
        #     COMPUTE_MEMORY_METRICS = [ 'Self GINTOP',
        #                                'Total GINTOP',
        #                                'Self GFLOP',
        #                                'Total GFLOP',
        #                                'Self Memory GB',
        #                                'Total Memory GB', ]
        #     loopFound = False
        #     with open(csvReport, 'r') as csvFile:
        #         # Skipping first 6 irrelevant lines..
        #         for i in range(6):
        #             next(csvFile)
        #         reader = csv.DictReader(csvFile)
        #         for row in reader:
        #             if row['Source Location'] == LOOP_SOURCE_LOCATION:
        #                 loopFound = True
        #                 for name in COMPUTE_MEMORY_METRICS:
        #                     if not name in row or not row[name]:
        #                         self.log.to_log((ICL_CPU_ROOFLINE_ANALYSIS + "FAILED: No or empty '{}'" \
        #                             + " for C++ loop").format(name))
        #                         self.status = 'FAIL'
        #                         status = 'Fail'
        #                         return title, status, important_messages
        #                     else:
        #                         self.log.to_log((ICL_CPU_ROOFLINE_ANALYSIS + ": Found non-empty '{}'" \
        #                             + " value for C++ loop: {}").format(name, row[name]))
        #                 break
        #     if not loopFound:
        #         self.log.to_log(ICL_CPU_ROOFLINE_ANALYSIS + "FAILED: Can't find C++ loop in the Survey report")
        #         self.status = 'FAIL'
        #         status = 'Fail'
        #         return title, status, important_messages

        status = 'Ok'
        return title, status, important_messages

    def get_errors_from_stderr(self, stderr, show_warnings=True):
        important_messages = []
        for line in stderr:
            if line.startswith(ERROR_TAG):
                important_messages.append(line)
            elif (show_warnings and line.startswith(WARNING_TAG)):
                # Not show unuseful warning for user
                if not CANNOT_LOCATE_FILE_WARNING in line and not CANNOT_LOCATE_DEBUGGING_INFORMATION_WARNING in line:
                    important_messages.append(line)
        
        return important_messages

    def remove_result_dir(self, test_descriptor):
        icl_cpu_result_dir_abs_path = self.state.append_to_work_dir(PROJECT_NAME)
        self.state.remove_dir(icl_cpu_result_dir_abs_path)

class AdvisorDpcppCPUTest(TestCase):
    def run(self):
        self.status = 'OK'
        self.log.to_stdout("=" * 80)
        self.log.to_stdout("Running DPC++ CPU test:")
        for test_descriptor in advisor_dpcpp_cpu_test_descriptor:
            self.log.to_log("=" * 80)
            test_title = test_descriptor['name']
            self.log.to_stdout(test_title)
            self.log.to_stdout(self.state.output_indent + 'Running collection...')

            title, status, important_messages = self.run_collection(test_descriptor)

            self.log.to_stdout(title + ': ' + status)

            for message in important_messages:
                self.log.to_stdout(message)

        self.log.to_stdout('')

        for test_descriptor in advisor_dpcpp_cpu_test_descriptor:
            if self.status == 'OK' and test_descriptor['report'] != '':
                report_title = test_descriptor['name']
                self.log.to_log("-" * 80)
                self.log.to_stdout(report_title + ' report')
                title, status, important_messages = self.run_report(test_descriptor)
                self.log.to_stdout(title + ': ' + status)

                for message in important_messages:
                    self.log.to_stdout(message)

        if self.status == 'OK':
            self.remove_result_dir(test_descriptor)
        
    def run_collection(self, test_descriptor):
        title = self.state.output_indent + COLLECTION

        args = [self.state.cl_path]
        args += ['--collect']
        args += test_descriptor['analysis']
        args += test_descriptor['analysis_opts']
        args += ['--project-dir']

        result_dir_abs_path = self.state.append_to_work_dir('dpcpp_cpu_result')

        args += [result_dir_abs_path]

        app_dir = os.path.join(self.state.product_dir, 'samples', 'en', 'C++', 'adv_check')
        app_name = DPCPP_CPU_EXE_NAME
        if sys.platform == 'win32':
            app_name += '.exe'
        app_path = os.path.join(app_dir, app_name)
        args += ['--', app_path]

        args += [DPCPP_CPU_EXE_OPTS]

        ret_code, stdout, stderr = self.subprocess_wrapper(args, indent=self.state.stdout_indent)

        important_messages = self.get_errors_from_stderr(stderr)

        if ret_code != 0:
            self.status = 'FAIL'
            status = 'Fail'
        else:
            status = 'Ok'

        return title, status, important_messages

    def run_report(self, test_descriptor):
        title = self.state.output_indent + REPORT

        result_dir_abs_path = self.state.append_to_work_dir('dpcpp_cpu_result')
        args = [self.state.cl_path,]
        args += ['--report']
        args += test_descriptor['report']
        args += ['--show-all-columns']
        args += ['--format', 'csv',]
        args += ['--project-dir']
        args += [result_dir_abs_path]

        ret_code, report, stderr = self.subprocess_wrapper(args, indent=self.state.stdout_indent)

        important_messages = self.get_errors_from_stderr(stderr)

        if ret_code != 0:
            self.status = 'FAIL'
            status = 'Fail'
            return title, status, important_messages

        tmp_str = ''.join(map(str, test_descriptor['report']))

        if test_descriptor['report'] == ['survey']:
                res_folder = 'hs000'

        csvReport = os.path.join(result_dir_abs_path, 'e000', res_folder, '{}-{}.csv'.format(CL_TOOL, tmp_str))

        if test_descriptor['report'] == ['survey']:
            LOOP_SOURCE_TO_CHECK = 'matrix.cpp:43'
            EXP_TRIP_COUNT = '512'
            EXP_CALL_COUNT = '262144'
            COMPUTE_MEMORY_METRICS = [ 'Self GINTOP',
                                       'Self GFLOP',
                                       'Self Memory GB', ]
            loopFound = False
            with open(csvReport, 'r') as csvFile:
                # Skipping first 6 irrelevant lines...
                for i in range(6):
                    next(csvFile)
                reader = csv.DictReader(csvFile)
                for row in reader:
                    if row['Source Location'] == LOOP_SOURCE_TO_CHECK:
                        loopFound = True
                        if row['Average Trip Count'] != EXP_TRIP_COUNT:
                            self.log.to_log(DPCPP_CPU_SURVEY_ANALYSIS + " FAILED: Incorrect \'Average Trip Count\' for DPC++ loop")
                            self.status = 'FAIL'
                            status = 'Fail'
                            return title, status, important_messages
                        else:
                            self.log.to_log(DPCPP_CPU_SURVEY_ANALYSIS + ": Found correct \'Average Trip Count\' for DPC++ loop: {}" \
                                .format(row['Average Trip Count']))
                        if row['Call Count'] != EXP_CALL_COUNT:
                            self.log.to_log(DPCPP_CPU_SURVEY_ANALYSIS + " FAILED: Incorrect \'Call Count\' for DPC++ loop")
                            self.status = 'FAIL'
                            status = 'Fail'
                            return title, status, important_messages
                        else:
                            self.log.to_log(DPCPP_CPU_SURVEY_ANALYSIS + ": Found correct \'Call Count\' for DPC++ loop: {}" \
                                .format(row['Call Count']))
                        for name in COMPUTE_MEMORY_METRICS:
                            if not name in row or not row[name]:
                                self.log.to_log((DPCPP_CPU_SURVEY_ANALYSIS + " FAILED: No or empty '{} for DPC++ loop'" \
                                    + " for C++ loop").format(name))
                                self.status = 'FAIL'
                                status = 'Fail'
                                return title, status, important_messages
                            else:
                                self.log.to_log((DPCPP_CPU_SURVEY_ANALYSIS + ": Found non-empty '{}' DPC++ loop" \
                                    + " value for C++ loop: {}").format(name, row[name]))
                        break
            if not loopFound:
                self.log.to_log(DPCPP_CPU_SURVEY_ANALYSIS + " FAIILED: Can't find DPC++ loop in the Survey report")
                self.status = 'FAIL'
                status = 'Fail'
                return title, status, important_messages

        status = 'Ok'
        return title, status, important_messages

    def get_errors_from_stderr(self, stderr, show_warnings=True):
        important_messages = []
        for line in stderr:
            if line.startswith(ERROR_TAG):
                important_messages.append(line)
            elif (show_warnings and line.startswith(WARNING_TAG)):
                # Not show unuseful warning for user
                if not CANNOT_LOCATE_FILE_WARNING in line and not CANNOT_LOCATE_DEBUGGING_INFORMATION_WARNING in line:
                    important_messages.append(line)
        
        return important_messages

    def remove_result_dir(self, test_descriptor):
        icl_cpu_result_dir_abs_path = self.state.append_to_work_dir('dpcpp_cpu_result')
        self.state.remove_dir(icl_cpu_result_dir_abs_path)


class AdvisorDpcppGPUTest(TestCase):
    def run(self):
        self.status = 'OK'
        self.log.to_stdout("=" * 80)
        self.log.to_stdout("Running DPC++ GPU test:")

        for test_descriptor in advisor_dpcpp_gpu_test_descriptor:
            if 'ONEAPI_DEVICE_SELECTOR' in os.environ:
                self.log.to_log('Unset ONEAPI_DEVICE_SELECTOR')
                del os.environ['ONEAPI_DEVICE_SELECTOR']
            if test_descriptor['name'] != DPCPP_GPU_L0_ROOFLINE_ANALYSIS:
                self.log.to_log('Set ONEAPI_DEVICE_SELECTOR to opencl')
                os.environ['ONEAPI_DEVICE_SELECTOR'] = 'opencl:gpu'

            test_title = test_descriptor['name']
            self.log.to_stdout(test_title)
            self.log.to_stdout(self.state.output_indent + 'Running collection...')

            title, status, important_messages = self.run_collection(test_descriptor)

            self.log.to_stdout(title + ': ' + status)

            for message in important_messages:
                self.log.to_stdout(message)

        self.log.to_stdout('')

        for test_descriptor in advisor_dpcpp_gpu_test_descriptor:
            if self.status == 'OK' and test_descriptor['report'] != '':
                report_title = test_descriptor['name']
                self.log.to_log("-" * 80)
                self.log.to_stdout(report_title + ' report')
                title, status, important_messages = self.run_report(test_descriptor)
                self.log.to_stdout(title + ': ' + status)

                for message in important_messages:
                    self.log.to_stdout(message)

        if self.status == 'OK':
            self.remove_result_dir(test_descriptor)

        if 'ONEAPI_DEVICE_SELECTOR' in os.environ:
            self.log.to_log('Unset ONEAPI_DEVICE_SELECTOR')
            del os.environ['ONEAPI_DEVICE_SELECTOR']


    def run_collection(self, test_descriptor):
        title = self.state.output_indent + COLLECTION

        args = [self.state.cl_path]
        args += ['--collect']
        args += test_descriptor['analysis']
        args += test_descriptor['analysis_opts']
        args += ['--project-dir']

        if test_descriptor['name'] == DPCPP_GPU_OCL_ROOFLINE_ANALYSIS:
            result_dir_abs_path = self.state.append_to_work_dir('dpcpp_gpu_result/ocl')
        else:
            result_dir_abs_path = self.state.append_to_work_dir('dpcpp_gpu_result/l0')

        args += [result_dir_abs_path]

        app_dir = os.path.join(self.state.product_dir, 'samples', 'en', 'C++', 'adv_check')
        app_name = DPCPP_GPU_EXE_NAME
        if sys.platform == 'win32':
            app_name += '.exe'
        app_path = os.path.join(app_dir, app_name)
        args += ['--', app_path]

        args += [DPCPP_GPU_EXE_OPTS]

        ret_code, stdout, stderr = self.subprocess_wrapper(args, indent=self.state.stdout_indent)

        important_messages = self.get_errors_from_stderr(stderr)

        if ret_code != 0:
            self.status = 'FAIL'
            status = 'Fail'
        else:
            status = 'Ok'

        return title, status, important_messages

    def run_report(self, test_descriptor):
        title = self.state.output_indent + REPORT
        adv_dir = os.path.join(self.state.product_dir, 'bin64')
        self.log.to_log("Dumping DPC++ GPU metrics...")

        if test_descriptor['name'] == DPCPP_GPU_OCL_ROOFLINE_ANALYSIS:
            profileDir = self.state.append_to_work_dir('dpcpp_gpu_result/ocl')
            gpuMetricsReport = os.path.join(profileDir, 'report_dpcpp_ocl_gpu_metrics.txt')
        else:
            profileDir = self.state.append_to_work_dir('dpcpp_gpu_result/l0')
            gpuMetricsReport = os.path.join(profileDir, 'report_dpcpp_l0_gpu_metrics.txt')            

        gpuDatasetCmd = 'advisor-python {} {}' .format(os.path.join(adv_dir, 'survey_gpu_dataset_test.py'), profileDir)

        with open(gpuMetricsReport, 'w') as fd:
            try:
                subprocess.check_call(gpuDatasetCmd.split(), stdout=fd)
            except:
                self.log.to_log("Failed to dump GPU metrics")
                self.status = 'FAIL'
                status = 'Fail'
                return title, status, ''

        with open(gpuMetricsReport, 'r') as fd:
            gpuScriptOutputLines = fd.readlines()
            GPU_TASK_TO_CHECK    = 'matrix_mult_kernel_parallel_for'
            EXP_GLOBAL_SIZE      = '2 x 2'
            GPU_COMPUTE_METRICS  = [ 'elapsed_time',
                                # Floating-point operations
                                'gpu_compute_performance_gflop',
                                'gpu_compute_performance_gflops',
                                'gpu_compute_performance_fp_ai',
                                # GTI traffic
                                'gpu_memory_data_transferred_gb_read',
                                'gpu_memory_data_transferred_gb_write',
                                # L3 traffic
                                'l3_shader_data_transferred_gb',
                                # CARM traffic
                                'carm_traffic_gb',
                                # Integer operations (element access operations are integer)
                                'gpu_compute_performance_gintop',
                                'gpu_compute_performance_gintops',
                                'gpu_compute_performance_int_ai', ]
            gpuTaskFound = False

            for line in gpuScriptOutputLines:
                gpuTask = json.loads(line.replace('\'', '"').replace('\\', '\\\\'))
                if gpuTask['computing_task'] == GPU_TASK_TO_CHECK:
                    gpuTaskFound = True
                    if gpuTask['work_size_global'] != EXP_GLOBAL_SIZE:
                        self.log.to_log("Incorrect \'work_size_global\' for DPC++ GPU task")
                        self.status = 'FAIL'
                        status = 'Fail'
                        return title, status, ''
                    else:
                        self.log.to_log("Found correct \'work_size_global\' for DPC++ GPU task: " + gpuTask['work_size_global'])
                    for name in GPU_COMPUTE_METRICS:
                        if not name in gpuTask or not gpuTask[name]:
                            self.log.to_log("No or empty '" + name + "' for DPC++ GPU task")
                            self.status = 'FAIL'
                            status = 'Fail'
                            return title, status, ''                            
                        else:
                            self.log.to_log("Found non-empty '" + name + "' value for DPC++ GPU task: '" + gpuTask[name] + "'")
                    break
            if not gpuTaskFound:
                self.log.to_log("Can't find DPC++ GPU task")
                self.status = 'FAIL'
                status = 'Fail'
                return title, status, ''

        status = 'Ok'

        return title, status, ''

    def get_errors_from_stderr(self, stderr, show_warnings=True):
        important_messages = []
        for line in stderr:
            if line.startswith(ERROR_TAG):
                important_messages.append(line)
            elif (show_warnings and line.startswith(WARNING_TAG)):
                # Not show unuseful warning for user
                if not CANNOT_LOCATE_FILE_WARNING in line and not CANNOT_LOCATE_DEBUGGING_INFORMATION_WARNING in line:
                    important_messages.append(line)
        
        return important_messages

    def remove_result_dir(self, test_descriptor):
        dpcpp_gpu_result_dir_abs_path = self.state.append_to_work_dir('dpcpp_gpu_result')
        self.state.remove_dir(dpcpp_gpu_result_dir_abs_path)

def main():
    parser = OptionParser(usage="Usage: %prog [options] arg")
    parser.add_option("--log-dir", dest="log_dir", help=HELP_LOG_DIR)
    (options, args) = parser.parse_args()

    status = 0

    work_dir = ''

    if options.log_dir:
        if not os.path.exists(options.log_dir):
            sys.stdout.write("%s\n" % FOLDER_NOT_AVAILABLE)
            return 1
        work_dir = options.log_dir

    # bin_dir = os.path.dirname(os.path.realpath(__file__))
    # bin_dir poinst to the directory where the advisor binaries are located
    bin_dir = "/opt/intel/oneapi/advisor/2024.0/bin64"

    state = State(bin_dir, work_dir, PRODUCT_ABBR, INTERNAL_PREFIX, CL_TOOL, COPYRIGHTS_ADV)

    if state.status == 'OK':
        state.log.to_log("Check of files: Ok")
    else:
        state.log.to_log("Check of files: Fail")
        return 1

    try:
        # Write context values to log
        context_values_test = ContextValuesTest(state)
        context_values_test.run()
        if context_values_test.status == 'FAIL':
            status = 1

        # Run Advisor icl CPU test
        test = AdvisorIclCPUTest(state)
        test.run()
        if test.status == 'FAIL':
            status = 1

        # # Run Advisor dpcpp CPU test
        # test = AdvisorDpcppCPUTest(state)
        # test.run()
        # if test.status == 'FAIL':
        #     status = 1

        # # Run Advisor dpcpp GPU test
        # test = AdvisorDpcppGPUTest(state)
        # test.run()
        # if test.status == 'FAIL':
        #     status = 1
    
    except Exception as e:
        state.log.to_log("Exception: %s" % e)
        state.log.to_stdout('\n' + UNKNOWN_FAIL)
        return 1

    if status == 0:
        state.log.to_stdout(SYSTEM_READY)
        if state.has_warnings:
            state.log.to_stdout(SEE_WARNINGS)
    else:
        state.log.to_stdout(SYSTEM_NOT_READY)

    return status


if __name__ == "__main__":
    ret_code = main()
    sys.exit(ret_code)
