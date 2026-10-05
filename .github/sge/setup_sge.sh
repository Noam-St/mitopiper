#!/usr/bin/env bash
# Install and start a one-machine Son of Grid Engine cluster on a GitHub Actions Ubuntu runner,
# with a queue "all.q" and the "shared" parallel environment mitopiper's job headers use.
set -euxo pipefail

host=$(hostname)
# Grid Engine wants the host name to resolve to a real (non-loopback) address.
ip=$(hostname -I | awk '{print $1}')
sudo sed -i "/[[:space:]]$host\$/d" /etc/hosts
echo "$ip $host" | sudo tee -a /etc/hosts

sudo debconf-set-selections <<SEL
gridengine-common shared/gridenginemaster string $host
gridengine-common shared/gridenginecell string default
gridengine-common shared/gridengineconfig boolean true
gridengine-client shared/gridenginemaster string $host
gridengine-client shared/gridenginecell string default
gridengine-client shared/gridengineconfig boolean true
postfix postfix/main_mailer_type select No configuration
SEL
sudo apt-get update -qq
# a huge open-files limit crashes Grid Engine's setup
ulimit -n 1024
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y -qq gridengine-master gridengine-exec gridengine-client
dpkg -l | grep gridengine

export SGE_ROOT=/var/lib/gridengine SGE_CELL=default
echo "SGE_ROOT=$SGE_ROOT" >> "${GITHUB_ENV:-/dev/null}"
echo "SGE_CELL=$SGE_CELL" >> "${GITHUB_ENV:-/dev/null}"
echo "$host" | sudo tee $SGE_ROOT/$SGE_CELL/common/act_qmaster
sudo service gridengine-master restart || sudo /etc/init.d/gridengine-master restart
sleep 5
pgrep -a sge_qmaster

sudo -u sgeadmin qconf -am "$USER"
qconf -ah "$host" || true
qconf -as "$host" || true

tmp=$(mktemp -d)
cat > "$tmp/host" <<HOST
hostname $host
load_scaling NONE
complex_values NONE
user_lists NONE
xuser_lists NONE
projects NONE
xprojects NONE
usage_scaling NONE
report_variables NONE
HOST
qconf -Ae "$tmp/host"
cat > "$tmp/allhosts" <<HG
group_name @allhosts
hostlist $host
HG
qconf -Ahgrp "$tmp/allhosts"
cat > "$tmp/pe" <<PE
pe_name            shared
slots              99
user_lists         NONE
xuser_lists        NONE
start_proc_args    NONE
stop_proc_args     NONE
allocation_rule    \$pe_slots
control_slaves     FALSE
job_is_first_task  TRUE
urgency_slots      min
accounting_summary FALSE
qsort_args         NONE
PE
qconf -Ap "$tmp/pe"
qconf -sq > /dev/null 2>&1 || true
# queue: start from the template, put every host in it, 4 slots, the "shared" PE
qconf -sq all.q > "$tmp/q" 2>/dev/null || sed -n '1,$p' /usr/share/gridengine/default/common/queue_template > "$tmp/q" 2>/dev/null || true
if [ ! -s "$tmp/q" ]; then
  cat > "$tmp/q" <<Q
qname all.q
hostlist @allhosts
seq_no 0
load_thresholds NONE
suspend_thresholds NONE
nsuspend 1
suspend_interval 00:05:00
priority 0
min_cpu_interval 00:05:00
processors UNDEFINED
qtype BATCH INTERACTIVE
ckpt_list NONE
pe_list shared
rerun FALSE
slots 4
tmpdir /tmp
shell /bin/bash
prolog NONE
epilog NONE
shell_start_mode posix_compliant
starter_method NONE
suspend_method NONE
resume_method NONE
terminate_method NONE
notify 00:00:60
owner_list NONE
user_lists NONE
xuser_lists NONE
subordinate_list NONE
complex_values NONE
projects NONE
xprojects NONE
calendar NONE
initial_state default
s_rt INFINITY
h_rt INFINITY
s_cpu INFINITY
h_cpu INFINITY
s_fsize INFINITY
h_fsize INFINITY
s_data INFINITY
h_data INFINITY
s_stack INFINITY
h_stack INFINITY
s_core INFINITY
h_core INFINITY
s_rss INFINITY
h_rss INFINITY
s_vmem INFINITY
h_vmem INFINITY
Q
  qconf -Aq "$tmp/q"
else
  qconf -mattr queue hostlist @allhosts all.q
  qconf -mattr queue slots 4 all.q
  qconf -mattr queue pe_list shared all.q
fi
sudo service gridengine-exec restart || sudo /etc/init.d/gridengine-exec restart
sleep 10
pgrep -a sge_execd
qhost
qstat -f
