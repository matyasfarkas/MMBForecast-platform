% Matlab/Octave toolbox for performing unit tests.
%
% Routines:
%
%   build_report_summary           - Builds a summary report.
%   dassert                        - Tests the equality of two objects.
%   dtest                          - Runs unit test defined in a Matlab/Octave routine, by calling mtest, and display results.
%   get_directory_description      - Lists recursively all the *.m files in a directory.
%   initialize_unit_tests_toolbox  - Initialization of the path to the m-unit-tests/src folder.
%   is_unitary_test_available      - Decides if unitary tests defined in a Matlab/Octave routine have to be run
%   mtest                          - Extracts unit test sections from Matlab/Octave's routine, executes the tests and reports results.
%   run_unitary_tests              - Runs unitary tests defined in a collection of files.
%   run_unitary_tests_in_directory - Runs all the unitary tests defined in a directory (and subfolders).
%
%
% Example:
%
% >> addpath m-unit-tests
% >> initialize_unit_test_toolbox
% >> run_unitary_tests_in_directory('dates/src/@dates')