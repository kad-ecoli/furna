#!/bin/bash

# This folder is the temporary folder for writing pdb files by cgi script.
# Permission and context of the folder and cgi scripts must be correctly set:

# First, set the permission:

    chmod a+x *.cgi
    chmod 777 output/

# Second, set the context, which may not be necessary on some system:

    chcon -t httpd_sys_script_exec_t *.cgi
    chcon -t httpd_sys_rw_content_t  output/
