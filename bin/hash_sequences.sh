#!/bin/bash

awk \
  '{ print $0 "\t" encodeData( $2 ) } 
  function encodeData( fld, cmd, output ) {
      cmd = "printf \047" fld "\047 | sha1sum"
      if ( (cmd | getline output) > 0 ) {
          sub(/ .*/,"",output)
      }
      else {
          print "failed to hash " fld | "cat>&2"
          output = fld
      }
      close( cmd )
      return output
  }' \
  "$1"
