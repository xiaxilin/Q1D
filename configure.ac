# -*- Autoconf -*-

AC_PREREQ(2.53)

define([git_describe], esyscmd([sh -c "git describe --all --tag | tr -d '\n'"]))dnl

AC_INIT([q1d], [git_describe], [joseph.derlaga@gmail.com])
AC_CONFIG_SRCDIR([q1d.f90])
AM_INIT_AUTOMAKE

AC_CONFIG_FILES([
  Makefile
])
AC_OUTPUT

