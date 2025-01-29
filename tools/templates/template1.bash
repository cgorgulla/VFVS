#!/usr/bin/env bash

# Copyright (C) 2019 Christoph Gorgulla
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
#
# This file is part of VirtualFlow.
#
# VirtualFlow is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# VirtualFlow is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with VirtualFlow.  If not, see <https://www.gnu.org/licenses/>.


# Job Information -- generally nothing in this
# section should be changed
##################################################################################

export VFVS_WORKUNIT={{workunit_id}}
export VFVS_JOB_STORAGE_MODE={{job_storage_mode}}
export VFVS_TMP_PATH=/dev/shm
export VFVS_CONFIG_JOB_TGZ={{job_tgz}}
export VFVS_TOOLS_PATH=${PWD}/bin
export VFVS_VCPUS={{threads_to_use}}

##################################################################################

for i in `seq 0 {{array_end}}`; do
	export VFVS_WORKUNIT_SUBJOB=${i}
	echo "Workunit ${VFVS_WORKUNIT}:${VFVS_WORKUNIT_SUBJOB}: stdout in {{batch_workunit_base}}/${VFVS_WORKUNIT_SUBJOB}.out, stderr in {{batch_workunit_base}}/${VFVS_WORKUNIT_SUBJOB}.err"
	date +%s > {{batch_workunit_base}}/${VFVS_WORKUNIT_SUBJOB}.start
	./templates/vfvs_run.py > {{batch_workunit_base}}/$$_${VFVS_WORKUNIT_SUBJOB}.out 2> {{batch_workunit_base}}/$$_${VFVS_WORKUNIT_SUBJOB}.err
	date +%s > {{batch_workunit_base}}/${VFVS_WORKUNIT_SUBJOB}.end
done
