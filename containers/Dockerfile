FROM mambaorg/micromamba:1.5.10
USER root
COPY environment.yml /opt/satpipe/environment.yml
RUN micromamba install -y -n base -f /opt/satpipe/environment.yml && micromamba clean -a -y
# the driver is bash+gawk (no pip package): the repository is copied to /opt/satpipe/src
COPY . /opt/satpipe/src
ENV PATH=/opt/satpipe/src/bin:/opt/conda/bin:$PATH LC_ALL=C
RUN gawk 'BEGIN{PROCINFO["sorted_in"]="@ind_str_asc"; print "gawk ok"}' && satpipe --version
ENTRYPOINT ["satpipe"]
