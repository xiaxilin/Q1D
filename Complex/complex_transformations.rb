require 'English'

module ComplexTransformations

  @@intrinsic_name_mangling = '%s_c'
  @@complex_functions_used = []

  INTRINSICS = %w( abs acos asin atan atan2
                   ceiling cosh cpu_time dim log10
                   max maxloc maxval min minloc minval nint
                   random_number sign tanh ).join('|')

  def ComplexTransformations.intrinsic_name_mangling=(new_mangling)
    @@intrinsic_name_mangling = new_mangling
  end

  def ComplexTransformations.complex_functions_used
    @@complex_functions_used.sort.uniq
  end

  def ComplexTransformations.complex_functions_used=(list)
    @@complex_functions_used = list
  end

  def full_transformation
    eval methods.grep(/.*complex.*/).sort.join('.')
  end

  def limited_transformation
    complexify_module_definition.complexify_module_usage
  end

  def insert_use_complex # stop. turn around. run.
    return self if @@complex_functions_used.empty?
    join.sub(/^ *(?!(!|'|"))(|\w+\s+)(program|module|subroutine|function)\s+\w+(\s*\([\w,\s\n&]*\))?/i) do |match|
      match << "\n  use complex_functions_c, only: " <<
        @@complex_functions_used.uniq.sort.join(
            ", &\n                                 ")
    end.split(/^/)
  end

  def complexify_real_variables
    prev_line_continued = false
    map do |line|
      unless prev_line_continued
        line.sub!(/^(\s*)real(?=\s*( |:|,|\((?!\s*system_)))/i,'\1complex' )
      end
      prev_line_continued = line.match(/[^!]*&/)
      line
    end
  end

  def complexify_module_definition
    map do |line|
      line.sub( /^(\s*)(end\s*|)(module\s+)(?!procedure)(\S+)/i, '\1\2\3\4_c' )
    end
  end

  def complexify_module_usage
    map do |line|
      line.sub( /^(\s*use\s+\w+)/i, '\1_c' )
    end
  end

  def complexify_include_file_name
    map do |line|
      line.sub( /^(\s*include\s+)(["'])(.+?)(\.f90)\2/i, '\1\2\3_c\4\2' )
    end
  end

  def find_complex_operators_used
    operators = inject([]) do |ops,line|
      ops << line.scan( /(>=|<=|>|<)/ )
    end.flatten.uniq.map{ |op| "operator(#{op})" }
    @@complex_functions_used = (@@complex_functions_used | operators).sort
    self
  end

  def complexify_intrinsics
    map do |line|
      line.gsub( /\b(#{INTRINSICS})(?=\s*\()/i ) do |match|
        intrinsic = @@intrinsic_name_mangling % $1
        if $PREMATCH.match(/!/) then
          match
        else
          @@complex_functions_used << intrinsic.downcase
          intrinsic
        end
      end
    end
  end

  def fun3d_activate_complex_mode
    map do |line|
      line.sub(/(\s*complex_mode\s*=\s*)(\.false\.)/i,'\1.true.')
    end
  end

  def fun3d_activate_flow_join
    map do |line|
      line.sub( /use io_c,\s*only\s*:\s*readme/,
          'use complex_helper, only : readme=>readme_driver' )
    end
  end

end
