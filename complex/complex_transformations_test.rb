#!/usr/bin/env ruby

require 'test/unit'
require 'complex_transformations'

class Array
  include ComplexTransformations
end

class TestInsertsUseComplexFunctions < Test::Unit::TestCase
  def setup
    ComplexTransformations.complex_functions_used = []
  end
  def test_do_not_insert_use_complex_functions_if_none_used
    assert_equal( ['module stuff'],
                  ['module stuff'].insert_use_complex )
    assert_equal( ['program stuff'],
                  ['program stuff'].insert_use_complex )
  end
  def test_insert_use_complex_functions_with_intrinsics
    assert_equal(
      ["module stuff\n","  use complex_functions_c, only: max_c\n"," max_c(1,2)"],
      ["module stuff\n max(1,2)"].complexify_intrinsics.insert_use_complex )
  end
  def test_insert_use_complex_functions_with_intrinsics_with_trick_module_name
    assert_equal(
      ["module a_module\n","  use complex_functions_c, only: min_c\n"," min_c(a,b)"],
      ["module a_module\n min(a,b)"].
      complexify_intrinsics.insert_use_complex )
  end
  def test_insert_complex_functions_if_subroutine
    assert_equal(
      ["subroutine s(args)\n","  use complex_functions_c, only: max_c\n"," max_c(1,2)"],
      ["subroutine s(args)\n max(1,2)"].complexify_intrinsics.insert_use_complex )
  end
  def test_insert_complex_functions_if_pure_function
    assert_equal(
      [" pure function f(a,b)\n","  use complex_functions_c, only: max_c\n"," max_c(1,2)"],
      [" pure function f(a,b)\n max(1,2)"].complexify_intrinsics.insert_use_complex )
  end
  def test_insert_complex_functions_with_continued_parameter_list
    assert_equal(
      ["subroutine s(a,&\n"," b)\n","  use complex_functions_c, only: max_c\n"," max_c(1,2)"],
      ["subroutine s(a,&\n b)\n max(1,2)"].complexify_intrinsics.insert_use_complex )
  end
end

class TestIntrinsicTranslation < Test::Unit::TestCase
  def setup
    ComplexTransformations.complex_functions_used = []
  end
  def test_translates_single_intrinsic_function
    assert_equal( ['a = abs_c(var)'],
                  ['a = abs(var)'].complexify_intrinsics )
  end
  def test_records_intrinsics_found
    ['a = min(var,2_dp) + abs(-1.0_dp)'].complexify_intrinsics
    assert_equal( ['abs_c','min_c'],
                  ComplexTransformations.complex_functions_used )
  end
  def test_translates_single_intrinsic_function_with_space
    assert_equal( ['a = acos_c (var)'],
                  ['a = acos (var)'].complexify_intrinsics )
  end
  def test_translates_nested_different_intrinsic_functions
    assert_equal( ['a = acos_c(min_c(var1,var2)))'],
                  ['a = acos(min(var1,var2)))'].complexify_intrinsics )
  end
  def test_translates_nested_same_intrinsic_functions
    assert_equal( ['a = min_c(min_c(var1,var2),var3))'],
                  ['a = min(min(var1,var2),var3))'].complexify_intrinsics )
  end
  def test_does_not_translates_words_containing_intrinsic_functions
    assert_equal( ['absolutely fabulous'],
                  ['absolutely fabulous'].complexify_intrinsics )
  end
  def test_ignore_variable_names_that_look_like_intrinsic_names
    assert_equal( ['x=ABS_c(sign_c(sign))'],
                  ['x=ABS(sign(sign))'].complexify_intrinsics )
  end
  def test_handles_parentheses_in_front_of_intrisic
    assert_equal( ['if ( real(abs_c(y(node1) - y(node2))) > y_coplanar_tol )'],
                  ['if ( real(abs(y(node1) - y(node2))) > y_coplanar_tol )'].
                  complexify_intrinsics )
  end
  def test_assume_that_real_is_not_an_intrinsic
    assert_equal( ['write(*,*) real(var)'],
                  ['write(*,*) real(var)'].complexify_intrinsics )
  end
  def test_should_not_convert_commented_intrinsic
    assert_equal( ['! c = min(a,b)'],
                  ['! c = min(a,b)'].complexify_intrinsics )
  end
  def test_does_not_translates_variable_names_containing_intrinsic_functions
    assert_equal( ['    rmax(:) = -huge(1.0_my_r8)'],
                  ['    rmax(:) = -huge(1.0_my_r8)'].complexify_intrinsics )
    assert_equal( ['    r_max(:) = -huge(1.0_my_r8)'],
                  ['    r_max(:) = -huge(1.0_my_r8)'].complexify_intrinsics )
    assert_equal( ['      if(rcheck(1) > rmax(1)) rmax(1) = rcheck(1)'],
                  ['      if(rcheck(1) > rmax(1)) rmax(1) = rcheck(1)'].
                  complexify_intrinsics )
  end
  def test_intrinsic_subroutine_call
    assert_equal( ['call cpu_time_c(starttime)'],
                  ['call cpu_time(starttime)'].complexify_intrinsics )
  end
  def test_intrinsic_subroutine_call_with_old_name_mangling
    begin
      ComplexTransformations.intrinsic_name_mangling = 'cc%s'
      assert_equal( ['call cccpu_time(starttime)'],
                    ['call cpu_time(starttime)'].complexify_intrinsics )
    ensure
      ComplexTransformations.intrinsic_name_mangling = '%s_c'
    end
  end
end

class TestFindComplexOperatorsUsed < Test::Unit::TestCase
  def setup
    ComplexTransformations.complex_functions_used = []
  end
  def test_finds_complex_operator_used
    ['a>1'].find_complex_operators_used
    assert_equal( [ 'operator(>)' ],
                  ComplexTransformations.complex_functions_used )
  end
  def test_finds_complex_operators_used
    ['b<0 .or. b>=1'].find_complex_operators_used
    assert_equal( [ 'operator(<)','operator(>=)' ],
                  ComplexTransformations.complex_functions_used )
  end
  def test_finds_lessthanequal_complex_operators_used
    ['b<0 .or. b<=1'].find_complex_operators_used
    assert_equal( [ 'operator(<)','operator(<=)' ],
                  ComplexTransformations.complex_functions_used )
  end
end

class TestRealVariableTranslations < Test::Unit::TestCase
  def test_translates_real_variables_to_complex
    assert_equal( ['complex variable'],
                  ['real variable'].complexify_real_variables )
    assert_equal( ['complex, dimension(:) :: variable'],
                  ['real, dimension(:) :: variable'].complexify_real_variables )
    assert_equal( ['complex::variable'],
                  ['real::variable'].complexify_real_variables )
    assert_equal( ['complex(dp) :: variable'],
                  ['real(dp) :: variable'].complexify_real_variables )
  end
  def test_preserves_indentation_when_translating_real_variables_to_complex
    assert_equal( [' complex variable'],
                  [' real variable'].complexify_real_variables )
  end
  def test_does_not_translates_real_in_comments_or_string
    assert_equal( ['complex(dp) :: var_x ! with real comment'],
                  ['real(dp) :: var_x ! with real comment'].
                  complexify_real_variables )
    assert_equal( ['!this is real cool'],
                  ['!this is real cool'].complexify_real_variables )
    assert_equal( ['write(*,*)"real"'],
                  ['write(*,*)"real"'].complexify_real_variables )
  end
  def test_does_not_translate_system_kinds
    assert_equal( ['real(system_r4) :: var_x'],
                  ['real(system_r4) :: var_x'].complexify_real_variables )
    assert_equal( ['real(system_r8) :: var_x'],
                  ['real(system_r8) :: var_x'].complexify_real_variables )
  end
  def test_does_not_get_fooled_by_system_dimensions
    assert_equal( ['complex(dqp), dimension(1_system_i1) :: var'],
                  ['real(dqp), dimension(1_system_i1) :: var'].
                  complexify_real_variables )
  end
  def test_ignore_real_intrinsics
    assert_equal( ['write(*,*) real(var)'],
                  ['write(*,*) real(var)'].complexify_real_variables )
  end
  def test_ignore_realz_varible_name
    assert_equal( ['realz, bc )'],
                  ['realz, bc )'].complexify_real_variables )
    assert_equal( ['  realz, bc )'],
                  ['  realz, bc )'].complexify_real_variables )
  end
  def test_continutation_lines_do_not_mask_a_real_instrinsic
    assert_equal( ['&',' real(var)'], 
                  ['&',' real(var)'].complexify_real_variables )
    assert_equal( ['& ! jerk',' real(var)'], 
                  ['& ! jerk',' real(var)'].complexify_real_variables )
  end
end

class TestFun3dSpecificTranslations < Test::Unit::TestCase
  def test_activates_complex_mode_boolean
    assert_equal( ['complex_mode = .true.'],
                  ['complex_mode = .FALSE.'].fun3d_activate_complex_mode )
    assert_equal( [' logical :: complex_mode = .true.'],
                  [' logical :: complex_mode = .false.'].
                    fun3d_activate_complex_mode )
  end
  def test_leaves_active_complex_mode_alone
    assert_equal( ['complex_mode = .true.'],
                  ['complex_mode = .true.'].fun3d_activate_complex_mode )
  end
  def test_leave_complex_mode_alone_in_conditional_test
    assert_equal( ['if ( complex_mode == .false. )'],
                  ['if ( complex_mode == .false. )'].
                  fun3d_activate_complex_mode )
  end
end

class TestUseModuleNameTranslations < Test::Unit::TestCase
  def test_translates_used_module_name
    assert_equal( [' use real_c'],
                  [' use real'].complexify_module_usage )
  end
  def test_translates_used_module_name_containing_numbers
    assert_equal( ['  use fun3d_c'],
                  ['  use fun3d'].complexify_module_usage )
  end
  def test_translates_used_module_name_with_only_qualifier
    assert_equal( [' use REAL_c, only : var'],
                  [' use REAL, only : var'].complexify_module_usage )
  end
  def test_translates_used_module_name_with_only_qualifier_space_before_comma
    assert_equal( [' use real_c , only : var'],
                  [' use real , only : var'].complexify_module_usage )
  end
end

class TestIncludeFileNameTranslations < Test::Unit::TestCase
  def test_translates_included_file_name
    assert_equal( [" include 'dog_c.f90'"],
                  [" include 'dog.f90'"].complexify_include_file_name )
  end
end

class TestModuleNameTranslations < Test::Unit::TestCase
  def test_does_not_translate_module_procedure
    assert_equal( ['   module procedure name'],
                  ['   module procedure name'].complexify_module_definition )
    assert_equal( ['   module Procedure name'],
                  ['   module Procedure name'].complexify_module_definition )
  end
  def test_translates_module_name
    assert_equal( ['module real_c'],
                  ['module real'].complexify_module_definition )
    assert_equal( ['end module real_c'],
                  ['end module real'].complexify_module_definition )
    assert_equal( ['endmodule real_c'],
                  ['endmodule real'].complexify_module_definition )
  end
  def test_translates_module_name_with_leading_spaces
    assert_equal( ['  module real_c'],
                  ['  module real'].complexify_module_definition )
    assert_equal( ['   end module real_c'],
                  ['   end module real'].complexify_module_definition )
  end
  def test_does_not_rename_module_name_if_name_missing
    assert_equal( ['module '], ['module '].complexify_module_definition )
  end
  def test_translates_module_names_preserving_spacing
    assert_equal( ['module  original_c'], 
                  ['module  original'].complexify_module_definition )
    assert_equal( [' end  module '],
                  [' end  module '].complexify_module_definition )
    assert_equal( [' endmodule '],
                  [' endmodule '].complexify_module_definition )
  end
  def test_translates_module_names_with_trailing_comments
    assert_equal( ['end module real_c ! comment'],
                  ['end module real ! comment'].complexify_module_definition )
  end
  def test_does_not_translate_module_keyword_in_comment
    assert_equal( ['! here is a module keyword in comment'],
                  ['! here is a module keyword in comment'].
                  complexify_module_definition )
  end
  def test_translates_module_names_capitailization_not_withstanding
    assert_equal( ['MODULE CAPITAL_c'],
                  ['MODULE CAPITAL'].complexify_module_definition )
  end
end

class TestJoinPointTranslations < Test::Unit::TestCase
  def test_flow_join_point
     assert_equal( ['use complex_helper, only : readme=>readme_driver'],
                   ['use io_c, only : readme'].fun3d_activate_flow_join)
  end
end
