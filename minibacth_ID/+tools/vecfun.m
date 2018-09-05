function [varargout] = vecfun(func, varargin)
    varargin = tools.cellfun(@(x) tools.matrix_vector(x), varargin);
    varargout = cell(nargout, 1);
    [varargout{:}] = arrayfun(func, varargin{:}, 'UniformOutput', false);
end


