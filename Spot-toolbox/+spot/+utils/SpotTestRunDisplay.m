classdef SpotTestRunDisplay < TestRunMonitor
%CommandWindowTestRunDisplay Print test suite execution results to Command Window.
%   CommandWindowTestRunDisplay is a subclass of TestRunMonitor.  If a
%   CommandWindowTestRunDisplay object is passed to the run method of a
%   TestComponent, such as a TestSuite or a TestCase, it will print information
%   to the Command Window as the test run proceeds.
%
%   CommandWindowTestRunDisplay methods:
%       testComponentStarted  - Update Command Window display
%       testComponentFinished - Update Command Window display
%       testCaseFailure       - Log test failure information
%       testCaseError         - Log test error information
%
%   CommandWindowTestRunDisplay properties:
%       TestCaseCount         - Number of test cases executed
%       Faults                - Struct array of test fault info
%
%   See also TestRunLogger, TestRunMonitor, TestSuite

%   Steven L. Eddins
%   Copyright 2008-2009 The MathWorks, Inc.
    
    properties (SetAccess = private)
        %TestCaseCount - Number of test cases executed
        TestCaseCount
        
        %Faults - Struct array of test fault info
        %   Faults is a struct array with these fields:
        %       Type      - either 'failure' or 'error'
        %       TestCase  - the TestCase object that suffered the fault
        %       Exception - the MException thrown when the fault occurred
        Faults = struct('Type', {}, 'TestCase', {}, 'Exception', {});
    end
    
    properties (SetAccess = private, GetAccess = private)
        %InitialTic - Out of tic at beginning of test run
        InitialTic
        
        %ComponentTic - Time for a particular test
        ComponentTic
        
        %InitialComponent First test component executed
        %   InitialComponent is set to the first test component executed in the
        %   test run.  This component is saved so that the end of the test run
        %   can be identified.
        InitialComponent = []        
    end
        
    
    methods
        
        function testComponentStarted(self, component)
            %testComponentStarted Update Command Window display
            %    If the InitialComponent property is not yet set, 
            %    obj.testComponentStarted(component) sets the property and calls
            %    obj.testRunStarted(component).
            
            if isempty(self.InitialComponent)
                self.InitialComponent = component;
                self.testRunStarted(component);
            end
            self.ComponentTic = tic;
        end    
            
        function testComponentFinished(self, component, did_pass)
            %testComponentFinished Update Command Window display
            %    If component is a TestCase object, then 
            %    obj.testComponentFinished(component, did_pass) prints pass/fail
            %    information to the Command Window.
            %
            %    If component is the InitialComponent, then
            %    obj.testRunFinished(did_pass) is called.
            
            if isa(component, 'TestCase')
                self.TestCaseCount = self.TestCaseCount + 1;
                fprintf(' %3i: %-40s\t\t',self.TestCaseCount,component.Name);
                if did_pass
                   fprintf('[pass]');
                else
                   fprintf('[fail]');
                end
                fprintf('%7.2f sec\n',toc(self.ComponentTic))
            end
            
            if isequal(component, self.InitialComponent)
                self.testRunFinished(did_pass);
            end
        end
               
        function testCaseFailure(self, test_case, failure_exception)
            %testCaseFailure Log test failure information
            %    obj.testCaseFailure(test_case, failure_exception) logs the test
            %    case failure information.
            
            self.logFault('failure', test_case, ...
                failure_exception);
        end
        
        function testCaseError(self, test_case, error_exception)
            %testCaseError Log test error information
            %    obj.testCaseError(test_case, error_exception) logs the test
            %    case error information.
            
            self.logFault('error', test_case, ...
                error_exception);
        end
        
    end
    
    methods (Access = private)
        function testRunStarted(self, component)
            %testRunStarted Update Command Window display
            %    obj.testRunStarted(component) displays information about the test
            %    run to the Command Window.
            
            self.InitialTic = tic;
            self.TestCaseCount = 0;
            num_cases = component.numTestCases();
            if num_cases == 1
                str = 'case';
            else
                str = 'cases';
            end
            fprintf('Starting test run with %d test %s.\n', ...
                num_cases, str);
        end
        
        function testRunFinished(self, did_pass)
            %testRunFinished Update Command Window display
            %    obj.testRunFinished(component) displays information about the test
            %    run results, including any test failures, to the Command Window.
            
            if did_pass
                result = 'PASSED';
            else
                result = 'FAILED';
            end
            
            fprintf('\n%s in %.3f seconds.\n', result, toc(self.InitialTic));
            
            self.displayFaults();
        end
        

        
        function logFault(self, type, test_case, exception)
            %logFault Log test fault information
            %    obj.logFault(type, test_case, exception) logs test fault
            %    information. type is either 'failure' or 'error'. test_case is a
            %    TestCase object.  exception is an MException object.
            
            self.Faults(end + 1).Type = type;
            self.Faults(end).TestCase = test_case;
            self.Faults(end).Exception = exception;
        end
        
        function displayFaults(self)
            %displayFaults Display test fault info to Command Window
            %    obj.displayFaults() displays a summary of each test failure and
            %    test error to the command window.
            for k = 1:numel(self.Faults)
                faultData = self.Faults(k);
                message = regexp(faultData.Exception.message,'\n','split');
                [tmp,file] = fileparts(faultData.TestCase.Location);
                if strcmp(faultData.Type, 'failure')
                    str = 'Failure';
                else
                    str = 'Error';
                end
                fprintf('\n%8s: %20s (in %20s): %s\n',...
                   str, file, faultData.TestCase.Name, message{1});
                displayStack(filterStack(faultData.Exception.stack));
                %fprintf('\n%s\n', faultData.Exception.message);
                %                keyboard
                fprintf('\n');
            end
        end
        
    end
    
end

function displayStack(stack)
%displayStack Display stack trace from MException instance
%   displayStack(stack) prints information about an exception stack to the
%   command window. 

for k = 1:numel(stack)
    filename = stack(k).file;
    linenumber = stack(k).line;
    href = sprintf('matlab: opentoline(''%s'',%d)', filename, linenumber);
    fprintf('    %s at <a href="%s">line %d</a>\n', filename, href, linenumber);
end
end

function new_stack = filterStack(stack)
%filterStack Remove unmeaningful stack trace calls
%    new_stack = filterStack(stack) removes from the input stack trace calls
%    that are framework functions and methods that are not likely to be
%    meaningful to the user.

% Testing stack traces follow this common pattern:
%
% 1. The first function call in the trace is often one of the assert functions
% in the framework directory.  This is useful to see.
%
% 2. The next function calls are in the user-written test functions/methods and
% the user-written code under test.  These calls are useful to see.
%
% 3. The final set of function calls are methods in the various framework
% classes.  There are usually several of these calls, which clutter up the 
% stack display without being that useful.
%
% The pattern above suggests the following stack filtering strategy: Once the
% stack trace has left the framework directory, do not follow the stack trace back
% into the framework directory.

mtest_directory = fileparts(which('runtests'));
last_keeper = numel(stack);
have_left_mtest_directory = false;
for k = 1:numel(stack)
    directory = fileparts(stack(k).file);
    if have_left_mtest_directory
        if strcmp(directory, mtest_directory)
            % Stack trace has reentered mtest directory.
            last_keeper = k - 1;
            break;
        end
    else
        if ~strcmp(directory, mtest_directory)
            have_left_mtest_directory = true;
        end
    end
end

new_stack = stack(1:last_keeper);
            
end

