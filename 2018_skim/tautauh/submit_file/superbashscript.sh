
usage() {
  echo -n "

This is my script template.

 Options:
  -u, --username    Username for script
  -p, --password    User password
  --force           Skip all user interaction.  Implied 'Yes' to all actions.
  -q, --quiet       Quiet (no output)
  -l, --log         Print log to file
  -s, --strict      Exit script with null variables.  i.e 'set -o nounset'
  -v, --verbose     Output more information. (Items echoed to 'verbose')
  -d, --debug       Runs script in BASH debug mode (set -x)
  -h, --help        Display this help and exit
      --version     Output version information and exit
"
}

# translate long options to short
for arg
do
    delim=""
    case "$arg" in
       --help) args="${args}-h ";;
       --verbose) args="${args}-v ";;
       --config) args="${args}-c ";;
       # pass through anything else
       *) [[ "${arg:0:1}" == "-" ]] || delim="\""
           args="${args}${delim}${arg}${delim} ";;
    esac
done
# reset the translated args
eval set -- $args
# now we can process with getopt
while getopts ":hvc:" opt; do
    case $opt in
        h)  usage ;;
        v)  echo "VERBOSE=true" ;;
        c)  echo "source $OPTARG" ;;
        \?) echo "usage" ;;
        :)
        echo "option -$OPTARG requires an argument"
        
        ;;
    esac
done
