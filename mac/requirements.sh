#!/bin/bash

set -xue

./src/init_lua-npge.sh

brew update
brew install homebrew/dupes/ncurses

for pkg in \
    qt4 cmake libzip libzzip \
    boost lua51 luabind \
    pandoc gnu-tar coreutils binutils \
; do
    # brew fails if a package is already installed:
    # https://travis-ci.org/npge/npge/jobs/151420671#L396
    # See also https://docs.travis-ci.com/user/osx-ci-environment/#A-note-on-upgrading-packages
    brew unlink $pkg || true
    brew install $pkg
done

ln -s /usr/local/bin/gtac /usr/local/bin/tac
