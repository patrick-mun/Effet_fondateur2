#!/usr/bin/env bash
set -euo pipefail

repository_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
conda_command="${CONDA_EXE:-/opt/anaconda3/condabin/conda}"
environment_name="effet-fondateur-ibd"
environment_prefix="/opt/anaconda3/envs/${environment_name}"
cache_root="${repository_root}/data/cache/tools/explicit_ibd"
hap_version="1.0.0-15Jun23.92f"
hap_commit="a3fdd2d6294903387ee29d2f67ce6260593b681f"
hap_sha256="52ba95a2a8990d212e53084ba2e017227e25d4357fb9a99f3cf182591946d34d"
refined_version="17Jan20.102"
refined_sha256="a2c2d7ee6c1dff5c06831a5987ba3a3ce831cb5798a13d413ab0c10a577db919"

if [[ ! -x "${environment_prefix}/bin/java" ]]; then
  "${conda_command}" create -y -n "${environment_name}" openjdk=17
fi

hap_dir="${cache_root}/hap-ibd/${hap_version}"
hap_source="${hap_dir}/source"
mkdir -p "${hap_dir}"
if [[ ! -d "${hap_source}/.git" ]]; then
  git clone https://github.com/browning-lab/hap-ibd.git "${hap_source}"
fi
git -C "${hap_source}" fetch --tags --force
git -C "${hap_source}" checkout --detach "${hap_commit}"
"${environment_prefix}/bin/javac" -cp "${hap_source}/src" "${hap_source}/src/hapibd/HapIbdMain.java"
"${environment_prefix}/bin/jar" --create --file "${hap_dir}/hap-ibd.jar" \
  --main-class hapibd.HapIbdMain --date=2023-06-15T00:00:00Z \
  -C "${hap_source}/src" ./

refined_dir="${cache_root}/refined-ibd/${refined_version}"
mkdir -p "${refined_dir}"
curl --fail --location --silent --show-error \
  "https://faculty.washington.edu/browning/refined-ibd/refined-ibd.${refined_version}.jar" \
  --output "${refined_dir}/refined-ibd.${refined_version}.jar"

printf '%s  %s\n' "${hap_sha256}" "${hap_dir}/hap-ibd.jar" | shasum -a 256 --check
printf '%s  %s\n' "${refined_sha256}" "${refined_dir}/refined-ibd.${refined_version}.jar" | shasum -a 256 --check
"${environment_prefix}/bin/java" -version
"${environment_prefix}/bin/java" -jar "${hap_dir}/hap-ibd.jar" 2>&1 | sed -n '1,2p'
"${environment_prefix}/bin/java" -jar "${refined_dir}/refined-ibd.${refined_version}.jar" 2>&1 | sed -n '1,2p'
