classdef MultivariateNormal < ContinuousDistribution
    % MultivariateNormal Multivariate Normal (Gaussian) Distribution
    %
    % Represents a d-dimensional normal distribution with mean vector mu and
    % covariance matrix Sigma. The distribution can be used standalone or
    % within a Prior for mixture models.
    %
    % Constructor:
    %   mvn = MultivariateNormal(mu, Sigma)
    %
    % Args:
    %   mu    - d-dimensional mean vector (column vector)
    %   Sigma - d x d positive definite covariance matrix
    %
    % Examples:
    %   % Create 2D normal distribution
    %   mu = [1; 2];
    %   Sigma = [1, 0.5; 0.5, 1];
    %   mvn = MultivariateNormal(mu, Sigma);
    %
    %   % Generate samples
    %   samples = mvn.sample(100);  % 100 x 2 matrix
    %
    %   % Evaluate PDF
    %   pdf_val = mvn.evalPDF([1; 2]);
    %
    %   % Extract marginal
    %   norm = mvn.getMarginalUniv(1);
    %
    % Copyright (c) 2012-2026, Imperial College London

    properties
        dimension;  % Dimensionality
    end

    methods
        function self = MultivariateNormal(mu, Sigma)
            % Validate inputs
            if ~isvector(mu)
                line_error(mfilename, 'mu must be a vector');
            end
            if ~ismatrix(Sigma) || size(Sigma,1) ~= size(Sigma,2)
                line_error(mfilename, 'Sigma must be square matrix');
            end

            % Force column vector
            mu = mu(:);
            d = length(mu);

            % Check dimension consistency
            if size(Sigma,1) ~= d
                line_error(mfilename, 'Sigma dimensions must match mu length');
            end

            % Check positive definite
            [~, p] = chol(Sigma);
            if p ~= 0
                line_error(mfilename, 'Sigma must be positive definite');
            end

            % Call superclass constructor
            self@ContinuousDistribution('MultivariateNormal', 2, [-Inf, Inf]);

            % Store dimension and parameters
            self.dimension = d;
            self.setParam(1, 'mu', mu);
            self.setParam(2, 'Sigma', Sigma);

            % Set mean of first component for Prior compatibility
            self.mean = mu(1);
        end

        % =================== ACCESSORS ===================

        function d = getDimension(self)
            % Get the dimensionality of the distribution
            d = self.dimension;
        end

        function mu = getMeanVector(self)
            % Get the mean vector (d x 1)
            mu = self.getParam(1).paramValue;
        end

        function Sigma = getCovariance(self)
            % Get the covariance matrix (d x d)
            Sigma = self.getParam(2).paramValue;
        end

        function R = getCorrelation(self)
            % Get the correlation matrix
            Sigma = self.getCovariance();
            d = self.dimension;

            % R(i,j) = Sigma(i,j) / (sqrt(Sigma(i,i)) * sqrt(Sigma(j,j)))
            R = zeros(d, d);
            for i = 1:d
                for j = 1:d
                    std_i = sqrt(Sigma(i,i));
                    std_j = sqrt(Sigma(j,j));
                    if std_i > eps && std_j > eps
                        R(i,j) = Sigma(i,j) / (std_i * std_j);
                    else
                        R(i,j) = double(i == j);
                    end
                end
            end
        end

        % =================== DISTRIBUTION METHODS ===================

        function MEAN = getMean(self)
            % Get the mean of the first component (for Prior compatibility)
            mu = self.getMeanVector();
            MEAN = mu(1);
        end

        function VAR = getVar(self)
            % Get the variance of the first component
            Sigma = self.getCovariance();
            VAR = Sigma(1,1);
        end

        function SCV = getSCV(self)
            % Get squared coefficient of variation (not meaningful for multivariate)
            SCV = NaN;
            line_warning(mfilename, 'SCV not defined for multivariate distributions');
        end

        function SKEW = getSkewness(self)
            % Get skewness (multivariate normal is symmetric)
            SKEW = 0;
        end

        % =================== SAMPLING ===================

        function X = sample(self, n)
            % Generate n samples from the multivariate normal distribution
            %
            % Returns n x d matrix of samples

            if nargin < 2
                n = 1;
            end

            mu = self.getMeanVector();
            Sigma = self.getCovariance();
            d = self.dimension;

            % Use mvnrnd if available (Statistics Toolbox)
            if exist('mvnrnd', 'file') == 2
                X = mvnrnd(mu', Sigma, n);
            else
                % Manual implementation using Cholesky decomposition
                % X = mu + L*Z where L = chol(Sigma, 'lower'), Z ~ N(0,I)
                L = chol(Sigma, 'lower');
                Z = randn(d, n);
                X = (mu + L*Z)';
            end
        end

        % =================== PDF EVALUATION ===================

        function p = evalPDF(self, x)
            % Evaluate the multivariate normal PDF
            %
            % Args:
            %   x - d x 1 column vector or d x n matrix of points
            %
            % Returns:
            %   p - PDF value(s)

            mu = self.getMeanVector();
            Sigma = self.getCovariance();
            d = self.dimension;

            % Handle both column vector and matrix input
            if isvector(x)
                x = x(:)';  % Convert to row vector
            end

            if exist('mvnpdf', 'file') == 2
                % Use Statistics Toolbox if available
                p = mvnpdf(x, mu', Sigma);
            else
                % Manual calculation: f(x) = (2π)^(-d/2) |Σ|^(-1/2) exp(-0.5(x-μ)'Σ^(-1)(x-μ))
                n = size(x, 1);
                p = zeros(n, 1);

                invSigma = inv(Sigma);
                detSigma = det(Sigma);
                normConst = 1 / sqrt((2*pi)^d * detSigma);

                for i = 1:n
                    diff = x(i,:)' - mu;
                    p(i) = normConst * exp(-0.5 * diff' * invSigma * diff);
                end
            end
        end

        function Ft = evalCDF(self, t)
            % CDF is not well-defined for multivariate distributions
            line_error(mfilename, 'CDF is not defined for multivariate distributions');
        end

        function L = evalLST(self, s)
            % LST is not well-defined for multivariate distributions
            line_error(mfilename, 'LST is not defined for multivariate distributions');
        end

        % =================== MARGINAL DISTRIBUTIONS ===================

        function mvn_marg = getMarginal(self, indices)
            % Extract a marginal distribution for a subset of dimensions
            %
            % Args:
            %   indices - vector of dimension indices to keep (1-based)
            %
            % Returns:
            %   mvn_marg - MultivariateNormal for the marginal

            mu = self.getMeanVector();
            Sigma = self.getCovariance();

            % Extract marginal mean and covariance
            mu_marg = mu(indices);
            Sigma_marg = Sigma(indices, indices);

            % Create new MultivariateNormal for marginal
            mvn_marg = MultivariateNormal(mu_marg, Sigma_marg);
        end

        function norm = getMarginalUniv(self, index)
            % Extract a univariate marginal distribution
            %
            % Args:
            %   index - dimension index (1-based)
            %
            % Returns:
            %   norm - 1D MultivariateNormal distribution for that dimension

            mu = self.getMeanVector();
            Sigma = self.getCovariance();

            mean_marg = mu(index);
            var_marg = Sigma(index, index);

            norm = MultivariateNormal(mean_marg, var_marg);
        end

        % =================== SERIALIZATION ===================

        function s = toString(self)
            % Convert to string representation
            s = sprintf('jline.MultivariateNormal(d=%d)', self.dimension);
        end
    end

    methods (Static)
        function mvn = fitMeanAndCovariance(mu, Sigma)
            % Create a multivariate normal from mean and covariance
            %
            % Args:
            %   mu - mean vector
            %   Sigma - covariance matrix
            %
            % Returns:
            %   mvn - MultivariateNormal distribution

            mvn = MultivariateNormal(mu, Sigma);
        end
    end
end
