#!/bin/bash

## Script to collapse homopolymer stretches to a single letter

bioawk \
  '{ gsub(/[A]+/,"A");gsub(/[C]+/,"C");gsub(/[T]+/,"T");gsub(/[G]+/,"G");gsub(/[N]+/,"N") }1' \
  "$1"
