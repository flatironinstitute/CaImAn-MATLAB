% Test suite for run_CNMF_patches
% Use 'run(TestRunCNMFPatches)' to run the whole suite.
classdef TestRunCNMFPatches < matlab.unittest.TestCase

    properties
        dataset   % generate dataset to test with
        datasize  % dataset size
        patches   % patches coordinates
    end

    properties (ClassSetupParameter)
        sizY = struct('data2d', [128, 128, 500], 'data3d', [128, 128, 7, 100]);
        memory_mapped = struct('memmapped', true, 'inRAM', false);
    end

    properties (MethodSetupParameter)
        patch_size = struct('default', {[]}, 'small', 25);
    end

    properties (TestParameter)
        cluster_pixels = struct('cluster_on', true, 'cluster_off', false);
    end

    methods (TestClassSetup)
        function createDataset(testCase, sizY, memory_mapped)
            rng(1)  % fix random generator seed for reproducibility
            Y = randn(sizY);
            if memory_mapped
                Yr = reshape(Y, [], sizY(end));
                F_dark = min(Yr(:));
                filename = [tempname, '.mat'];
                save(filename, 'Yr', 'Y', 'F_dark', 'sizY', '-v7.3');
                Y = matfile(filename, 'Writable', false);
            end
            testCase.datasize = sizY;
            testCase.dataset = Y;
        end
    end

    methods (TestClassTeardown)
        function cleanDataset(testCase)
            % clear create file, if any
            if isa(testCase.dataset, 'matlab.io.MatFile')
                delete(testCase.dataset.Properties.Source);
            end
        end
    end

    methods (TestMethodSetup)
        function createPatches(testCase, patch_size)
            if isempty(patch_size)
                testCase.patches = [];
            else
                data_ndims = length(testCase.datasize(1:end-1));
                psize = [patch_size, patch_size, 5];
                psize = psize(1:data_ndims);
                testCase.patches = construct_patches( ...
                    testCase.datasize(1:end-1), psize);
            end
        end
    end

    methods (Test)
        function testRunning(testCase, cluster_pixels)
            % simply test that nothing crashes
            options = struct('cluster_pixels', cluster_pixels);
            run_CNMF_patches( ...
                testCase.dataset, 7, testCase.patches, [], [], options);
        end
    end
end