# Falcon C Utilities

## Requirements
We need **htslib** at runtime (via `dlopen()`).

All the **nim** dependencies are under `vendor/` as
"subtrees". (See `vendor/readme.md`)

## Set up

    source module.sh
    export NIMBLE_DIR=$(pwd)/.git/nimble  # recommended

## Test
To run unit-tests in the "tests/" directory:

    nimble test

## Test integration

    nimble integ  # redundant with our unit-tests, but you get the idea

## Install

    nimble install

## Debug and develop

    make
    make integ
