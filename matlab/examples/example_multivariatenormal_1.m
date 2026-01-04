function example_multivariatenormal_1()
    % Example demonstrating MultivariateNormal distribution usage
    %
    % This example shows:
    % 1. Creating a multivariate normal distribution
    % 2. Accessing statistical properties
    % 3. Sampling from the distribution
    % 4. Evaluating the PDF
    % 5. Extracting marginal distributions
    % 6. Using MultivariateNormal within a Prior for mixture models

    fprintf('=== MultivariateNormal Distribution Example ===\n\n');

    % Example 1: Basic 2D distribution
    example1_basic_2d_distribution();

    % Example 2: Sampling and statistics
    example2_sampling_and_statistics();

    % Example 3: PDF evaluation
    example3_pdf_evaluation();

    % Example 4: Marginal distributions
    example4_marginal_distributions();

    % Example 5: Mixture model with Prior
    example5_mixture_model_with_prior();

    fprintf('\n=== All examples completed successfully ===\n');
end

% =================== EXAMPLE 1: BASIC 2D DISTRIBUTION ===================

function example1_basic_2d_distribution()
    fprintf('Example 1: Basic 2D Distribution\n');
    fprintf('--------------------------------\n');

    % Create mean vector
    mu = [1; 2];

    % Create covariance matrix
    Sigma = [1.0, 0.5; 0.5, 1.0];

    % Create multivariate normal
    mvn = MultivariateNormal(mu, Sigma);

    % Inspect properties
    fprintf('Dimensionality: %d\n', mvn.getDimension());
    fprintf('Mean: [%.2f, %.2f]\n', mvn.getMeanVector()(1), mvn.getMeanVector()(2));
    fprintf('Variance (dim 1): %.4f\n', mvn.getVar());
    fprintf('Skewness: %.4f (always 0 for MVN)\n\n', mvn.getSkewness());
end

% =================== EXAMPLE 2: SAMPLING AND STATISTICS ===================

function example2_sampling_and_statistics()
    fprintf('Example 2: Sampling and Statistics\n');
    fprintf('-----------------------------------\n');

    mu = [0; 0];
    Sigma = [2.0, 0.8; 0.8, 1.0];
    mvn = MultivariateNormal(mu, Sigma);

    % Generate large sample
    n = 10000;
    samples = mvn.sample(n);

    % Compute sample statistics
    emp_mean = mean(samples);
    emp_cov = cov(samples);

    fprintf('True mean: [%.4f, %.4f]\n', mu(1), mu(2));
    fprintf('Sample mean: [%.4f, %.4f]\n\n', emp_mean(1), emp_mean(2));

    fprintf('True covariance:\n');
    fprintf('  [%.4f  %.4f]\n', Sigma(1,1), Sigma(1,2));
    fprintf('  [%.4f  %.4f]\n\n', Sigma(2,1), Sigma(2,2));

    fprintf('Sample covariance:\n');
    fprintf('  [%.4f  %.4f]\n', emp_cov(1,1), emp_cov(1,2));
    fprintf('  [%.4f  %.4f]\n\n', emp_cov(2,1), emp_cov(2,2));
end

% =================== EXAMPLE 3: PDF EVALUATION ===================

function example3_pdf_evaluation()
    fprintf('Example 3: PDF Evaluation\n');
    fprintf('-------------------------\n');

    % Standard 2D normal
    mvn = MultivariateNormal([0; 0], eye(2));

    % At mean: f(0,0) = 1/(2π) ≈ 0.1592
    pdf_at_mean = mvn.evalPDF([0; 0]);
    fprintf('PDF at mean (0,0): %.6f (expected: %.6f)\n', pdf_at_mean, 1/(2*pi));

    % Away from mean
    pdf_away = mvn.evalPDF([1; 1]);
    fprintf('PDF at (1,1): %.6f\n', pdf_away);

    % Far from mean
    pdf_far = mvn.evalPDF([3; 3]);
    fprintf('PDF at (3,3): %.6f\n\n', pdf_far);
end

% =================== EXAMPLE 4: MARGINAL DISTRIBUTIONS ===================

function example4_marginal_distributions()
    fprintf('Example 4: Marginal Distributions\n');
    fprintf('----------------------------------\n');

    % Create 3D distribution
    mu = [1; 2; 3];
    Sigma = [1.0, 0.5, 0.2; 0.5, 2.0, 0.3; 0.2, 0.3, 1.5];
    mvn3d = MultivariateNormal(mu, Sigma);

    fprintf('3D Distribution: dimension = %d\n\n', mvn3d.getDimension());

    % Extract 2D marginal (dimensions 1 and 3, using 1-based indexing)
    mvn2d = mvn3d.getMarginal([1, 3]);
    fprintf('2D Marginal (dims 1,3): dimension = %d\n', mvn2d.getDimension());
    mu_marg = mvn2d.getMeanVector();
    fprintf('Marginal mean: [%.2f, %.2f]\n\n', mu_marg(1), mu_marg(2));

    % Extract univariate marginals
    norm1 = mvn3d.getMarginalUniv(1);
    norm2 = mvn3d.getMarginalUniv(2);
    norm3 = mvn3d.getMarginalUniv(3);

    fprintf('Univariate marginals:\n');
    fprintf('  Dim 1: N(%.2f, %.4f)\n', norm1.getMean(), norm1.getVar());
    fprintf('  Dim 2: N(%.2f, %.4f)\n', norm2.getMean(), norm2.getVar());
    fprintf('  Dim 3: N(%.2f, %.4f)\n\n', norm3.getMean(), norm3.getVar());
end

% =================== EXAMPLE 5: MIXTURE MODEL WITH PRIOR ===================

function example5_mixture_model_with_prior()
    fprintf('Example 5: Mixture Model with Prior\n');
    fprintf('------------------------------------\n');

    % Create two alternative MVN distributions
    mvn1 = MultivariateNormal([0; 0], [1.0, 0.3; 0.3, 1.0]);
    mvn2 = MultivariateNormal([2; 2], [0.5, 0.1; 0.1, 0.5]);

    % Create Prior with 60% probability for mvn1, 40% for mvn2
    prior = Prior({mvn1, mvn2}, [0.6, 0.4]);

    fprintf('Prior mixture model:\n');
    fprintf('  Alternative 1 (60%%): MVN(μ=[0,0]'', Σ=[1,0.3;0.3,1])\n');
    fprintf('  Alternative 2 (40%%): MVN(μ=[2,2]'', Σ=[0.5,0.1;0.1,0.5])\n\n');

    % Properties of the mixture
    fprintf('Mixture properties:\n');
    fprintf('  Number of alternatives: %d\n', prior.getNumAlternatives());
    fprintf('  Mean (first component): %.4f\n', prior.getMean());
    fprintf('  (Expected: 0.6*0 + 0.4*2 = 0.8)\n\n');

    % Generate mixture samples
    fprintf('Generating 10 samples from the mixture:\n');
    samples = prior.sample(10);
    for i = 1:10
        fprintf('  Sample %d: %.4f\n', i, samples(i));
    end
end
