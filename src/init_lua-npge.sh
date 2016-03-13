#!/bin/bash

set -xue

# fill src/lua-npge with files from https://github.com/npge/lua-npge

COMMIT="388bae6600f47b0806797f90b56e8ae268463645"

lua_npge_dir="src/lua-npge"
lua_npge_git="$lua_npge_dir/.git"
git_cmd="git --git-dir=$lua_npge_git --work-tree=$lua_npge_dir"

if [ -e "$lua_npge_git" ]; then
    git submodule update --init
else
    # clone lua-npge from GitHub
    git clone https://github.com/npge/lua-npge $lua_npge_dir
    $git_cmd reset --hard $COMMIT
fi

# compare commits
actual_commit="$($git_cmd log --pretty=tformat:%H -1)"
if [ "$actual_commit" != "$COMMIT" ]; then
    echo "Bad lua-npge revision: $actual_commit (expected $COMMIT)"
    exit 1
else
    echo "OK. Good lua-npge revision"
fi
