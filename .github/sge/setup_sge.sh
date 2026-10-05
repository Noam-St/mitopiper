#!/usr/bin/env bash
# Build and start a one-machine Grid Engine cluster on an Ubuntu GitHub Actions runner:
# queue "all.q" on this host and the "shared" parallel environment mitopiper's job headers use.
#
# Ubuntu's gridengine packages (8.1.9, 22.04 and 24.04) crash while initialising the cluster
# (segfault in the spooling library), so Grid Engine is built from daimh/sge, a maintained fork
# of Son of Grid Engine that is tested weekly on Ubuntu. These are the steps of its README's
# "Quick test on one machine", checked by hand before being put here.
set -euxo pipefail

SGE_COMMIT=${SGE_COMMIT:-82bf084f03d3a12d4a2d55aea83ca7b1915e1c29}
host=$(hostname)
ip=$(hostname -I | awk '{print $1}')
sudo sed -i "/[[:space:]]$host\$/d" /etc/hosts      # the host name must not resolve to loopback
echo "$ip $host" | sudo tee -a /etc/hosts

sudo apt-get update -qq
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y -qq build-essential cmake git libdb5.3-dev \
    libhwloc-dev libncurses-dev libpam0g-dev libssl-dev libsystemd-dev libtirpc-dev libxext-dev pkgconf
src=$(mktemp -d)
git clone -q https://github.com/daimh/sge "$src"
git -C "$src" checkout -q "$SGE_COMMIT"
cmake -S "$src" -B "$src/build" -DCMAKE_INSTALL_PREFIX=/opt/sge > /dev/null
cmake --build "$src/build" -j "$(nproc)" > /dev/null
sudo cmake --install "$src/build" > /dev/null

id sge > /dev/null 2>&1 || sudo useradd -r -d /opt/sge sge
sudo chown -R sge /opt/sge
cd /opt/sge
# The installers also try to start the daemons through systemd; start them by hand when that fails.
export SGE_ROOT=/opt/sge SGE_CELL=default
(yes "" | sudo -E ./install_qmaster > /tmp/install_qmaster.log 2>&1) || true
tail -n 25 /tmp/install_qmaster.log
test -f /opt/sge/default/common/settings.sh      # the installer must at least have written the configuration
pgrep sge_qmaster > /dev/null || sudo -E /opt/sge/bin/lx-amd64/sge_qmaster
sleep 3
(yes "" | sudo -E ./install_execd > /tmp/install_execd.log 2>&1) || true
tail -n 15 /tmp/install_execd.log
pgrep sge_execd > /dev/null || sudo -E /opt/sge/bin/lx-amd64/sge_execd
sleep 5
pgrep -a sge_

set +u; source /opt/sge/default/common/settings.sh; set -u
{ echo "SGE_ROOT=$SGE_ROOT"; echo "SGE_CELL=$SGE_CELL"; } >> "${GITHUB_ENV:-/dev/null}"
echo "/opt/sge/bin/lx-amd64" >> "${GITHUB_PATH:-/dev/null}"

tmp=$(mktemp -d)
printf 'group_name @allhosts\nhostlist %s\n' "$host" > "$tmp/allhosts"
sudo -E env PATH="$PATH" qconf -Ahgrp "$tmp/allhosts"
qconf -sp smp | sed 's/^pe_name .*/pe_name            shared/' > "$tmp/pe"
sudo -E env PATH="$PATH" qconf -Ap "$tmp/pe"
cat > "$tmp/edit_queue" <<'EOF'
#!/bin/sh
sed -i -e 's/^qname .*/qname all.q/' -e 's/^hostlist .*/hostlist @allhosts/' -e 's/^slots .*/slots 4/' \
       -e 's/^pe_list .*/pe_list make smp shared/' -e 's/^shell .*/shell \/bin\/bash/' "$1"
EOF
chmod +x "$tmp/edit_queue"
sudo -E env PATH="$PATH" EDITOR="$tmp/edit_queue" qconf -aq all.q
sudo -E env PATH="$PATH" qconf -as "$host"
qhost -q
qstat -f
