#!/bin/bash

bioawk \
  '{ gsub(/[A]+/,"A");gsub(/[C]+/,"C");gsub(/[T]+/,"T");gsub(/[G]+/,"G");gsub(/[N]+/,"N") }1' \
  "$1"
