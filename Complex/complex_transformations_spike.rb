#! /usr/bin/env ruby
#
# Orchestrates forward-mode complexification on a file and
# directory level.
#
# This includes weaving complex aspects in at join points as
# appropriate to handle reading real or complex input depending
# on the situation.

$:.push File.dirname(__FILE__)

require 'complex_transformations'

LIMITED_TRANSFORMS = %w[
  complex_functions.f90
  set_precision.f90
]

class Array
  include ComplexTransformations
  @@intrinsic_name_mangling = 'cc%s'
end

files = ARGV.empty? ? Dir['*/*.?90'] : ARGV

files.delete_if{ |name| name.match(/_c\.(f|F)90$/i) }

files.each do |file|
  file_c = file.sub(/^(.*)\.(F|f)90$/,'\1_c.\290')
  next if File.exists?(file_c) && File.mtime(file_c) > File.mtime(file)
  ComplexTransformations.complex_functions_used = []
  lines = File.readlines(file)
  File.open(file_c,'w') do |f|
    if LIMITED_TRANSFORMS.include?(File.basename(file)) then
      puts "#{file_c}(limited transformation)"
      f.puts lines.limited_transformation
    else
      lines = lines.full_transformation
      path, filename = File.split(file)
      dir_file = File.join(File.basename(path),filename)
      case dir_file
      when 'FUN3D_90/main.F90'
        puts "#{file_c}(flow join point)"
        lines = lines.fun3d_activate_flow_join
      else
        puts file_c
      end
      f.puts lines
    end
  end
end
