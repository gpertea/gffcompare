#!/usr/bin/env bash

err_exit() {
  local msg="$1"
  local code=1 # Default exit code
  if [[ "$2" =~ ^[0-9]+$ ]]; then
    code="$2" # use it as the exit code.
  elif [[ -n "$2" ]]; then ## usage text?
    # This uses Bash's indirect parameter expansion. Output goes to stdout.
    printf '%s\n' "${!2}"
  fi
  # Print the final error message to stderr.
  echo -e "Error: $msg" >&2
  exit "$code"
}

if [[ ! -d tests/in ]]; then
 err_exit " tests/in not found!"
fi

mkdir -p tests/wrk

cd tests/wrk

## test 1 : self-matching check

../../gffcompare -T ../in/self-match.gtf -r ../in/self-match.gtf -o selfm

fexp=../expected_out/selfm.stats
fout=selfm.stats
if diff -q -I '^#' $fout $fexp &>/dev/null; then
    echo "  OK."
else
   err_exit "Error: test failed, output $fout different than expected ($fexp)!"
fi


