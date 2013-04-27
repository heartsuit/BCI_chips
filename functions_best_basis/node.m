function index = node(d,b)
% node -- Tree indexing function
%  Usage
%    index = node(d,b)
%  Inputs
%    d        depth from root of tree
%    b        index among the 2^d possibilities
%             in a left-right scan at that depth
%  Outputs
%    index    linear index of node in tree structure
%
	index =  2^d + b;
