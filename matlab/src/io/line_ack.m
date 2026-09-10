function line_ack(toolName, verbose, msg, cite)
% LINE_ACK(TOOLNAME, VERBOSE, MSG, CITE)
%
% Print, once per session, the acknowledgement of the external tool that a
% wrapper solver delegates to, together with the pointer to its official
% website and the canonical paper to cite. The acknowledgement is pull-based,
% like the library attribution: nothing is printed at the default verbosity,
% and only a caller that asks for it explicitly, by running at
% VerboseLevel.DEBUG, gets the line. It is printed at most once per TOOLNAME
% per session. SOLVER.CITATIONS and LINE_CITATION are the quiet ways to obtain
% the same reference.
%
% MSG and CITE supply the acknowledgement text and the reference for a solver
% that lives outside this tree; when omitted both come from the table below,
% which covers the in-tree wrapper solvers only. Mirror any edit to that table
% in the Java (jline.io.InputOutput.line_ack) and Python
% (line_solver.api.io.logging.line_ack) tables.
%
% The machine-readable form of the same reference is LINE_CITATION(TOOLNAME),
% which returns the BibTeX entry.
%
% See also LINE_CITATION.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    verbose = GlobalConstants.Verbose;
end

if verbose ~= VerboseLevel.DEBUG || GlobalConstants.Verbose == VerboseLevel.SILENT
    return
end

if GlobalConstants.isToolAckShown(toolName)
    return
end

if nargin < 3 || isempty(msg)
    msg = ackText(toolName);
end
if nargin < 4
    cite = '';
end
if isempty(cite)
    cite = ackCitation(toolName);
end
if isempty(msg)
    return
end

line_printf('%s\n', msg);
if ~isempty(cite)
    line_printf('  Cite: %s\n', cite);
end
GlobalConstants.setToolAckShown(toolName);
end

function msg = ackText(toolName)
% Attribution strings verified against each upstream project's own pages.
% Do not reword an author list or a URL from memory.
switch upper(toolName)
    case 'JMT'
        msg = ['SolverJMT delegates to Java Modelling Tools (JMT), by M. Bertoli, ', ...
            'G. Casale, G. Serazzi (Politecnico di Milano, Imperial College London). ', ...
            'Please acknowledge the JMT authors: http://jmt.sourceforge.net/'];
    case 'LQNS'
        msg = ['SolverLQNS delegates to LQNS/LQSIM, by G. Franks, M. Woodside et al. ', ...
            '(Real-Time and Distributed Systems Group, Carleton University). ', ...
            'Please acknowledge the LQNS authors: http://www.layeredqueues.org/'];
    case 'QNS'
        msg = ['SolverQNS delegates to qnsolver, part of the LQNS distribution by ', ...
            'G. Franks, M. Woodside et al. (Real-Time and Distributed Systems Group, ', ...
            'Carleton University). Please acknowledge the LQNS authors: ', ...
            'http://www.layeredqueues.org/'];
    otherwise
        % Out-of-tree wrapper solvers supply their own text through MSG: their
        % tools must not be named in this codebase.
        msg = '';
end
end

function cite = ackCitation(toolName)
% One-line reference to the canonical paper of each tool. It describes the same
% work as the BibTeX entry returned by LINE_CITATION (keys BerCS07 and
% FraAWDD09 in doc/latex/biblio.bib), which is where the citation key belongs:
% the printed line is for the reader, not for a .bib file.
switch upper(toolName)
    case 'JMT'
        cite = ['M. Bertoli, G. Casale, G. Serazzi. "The JMT Simulator for ', ...
            'Performance Evaluation of Non-Product-Form Queueing Networks". ', ...
            'Proc. of the 40th Annual Simulation Symposium (ANSS), pp. 3-10, 2007.'];
    case {'LQNS','QNS'}
        cite = ['G. Franks, T. Al-Omari, M. Woodside, O. Das, S. Derisavi. ', ...
            '"Enhanced Modeling and Solution of Layered Queueing Networks". ', ...
            'IEEE Trans. Software Engineering, 35(2):148-161, 2009.'];
    otherwise
        cite = '';
end
end
