tests_dir   := unit_tests
test_srcs   := $(wildcard $(source_dir)/$(tests_dir)/*.cpp)
test_objs   := $(patsubst %.cpp, %.o, $(subst $(source_dir),$(build_dir),$(test_srcs)))
test_target := $(build_dir)/$(tests_dir)/tests

test_depend_dir   := $(depend_dir)/$(tests_dir)
test_depend_files := $(wildcard $(test_depend_dir)/*$(depend_suffix))
unit-tests: depend_dir   := $(test_depend_dir)
unit-tests: depend_files := $(test_depend_files)

# Use a conditional include so consecutive cleans work
-include $(depend_files)

$(build_dir)/$(tests_dir)/%.o: $(source_dir)/$(tests_dir)/%.cpp
	$(compile-cxx)

unit-tests: $(test_target)

$(test_target): ld_libs := $(cuda_libs) $(structures_path)/libTipsy.a $(threads)
$(test_target): $(test_objs)
	$(link-charmc)

.PHONY: tests-clean
tests-clean:
	@ $(RM) $(filter-out $(build_dir)/$(tests_dir)/driver.o,$(test_objs))

.PHONY: tests-dist-clean
tests-dist-clean: tests-clean 
	@ $(RM) $(build_dir)/driver.o $(test_target) $(test_depend_files)

dist-clean: tests-dist-clean
clean: tests-clean
