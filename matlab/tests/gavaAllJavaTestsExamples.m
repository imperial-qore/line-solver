clc;
cwd = fileparts(mfilename('fullpath'))
cd(cwd);
outDir = fullfile(cwd,'generated');
if ~exist(outDir,'dir'), mkdir(outDir); end
fname = fullfile(outDir,'ExampleClosedModel.java');
fid = fopen(fname,'w+');

fprintf(fid, 'package jline.solvers;\n\n');
fprintf(fid, 'import jline.lang.*;\n');
fprintf(fid, 'import jline.lang.constant.*;\n');
fprintf(fid, 'import jline.lang.nodes.*;\n');
fprintf(fid, 'import jline.lang.distributions.*;\n');
fprintf(fid, 'import jline.lang.nodes.Queue;;\n');
fprintf(fid, 'import jline.solvers.mva.SolverMVA;\n\n');
fprintf(fid, 'import org.junit.jupiter.api.Test;\n\n');
fprintf(fid, 'import java.util.*;\n\n');
fprintf(fid, 'import static jline.examples.Examples.*;\n');
fprintf(fid, 'import static org.junit.jupiter.api.Assertions.assertEquals;\n\n');
fprintf(fid, 'public class ExampleClosedModel {\n\n');
fclose(fid);
%%
GlobalConstants.setDummyMode(true); exampleName='cqn_repairmen'; eval(exampleName); GlobalConstants.setDummyMode(false);
fname = fullfile(outDir,'ExampleClosedModel.java');
fid = fopen(fname,'a');
fprintf(fid,'\n\t@Test\n');
fprintf(fid,'\tpublic void test_%s() throws IllegalAccessException {\n',exampleName);
QN2JAVA(model, exampleName, fid, false);
printExample(fid, exampleName, model, solver, AvgTable);
fprintf(fid, '\t}\n');
fclose(fid);
%%
GlobalConstants.setDummyMode(true); exampleName='cqn_twoclass_hyperl'; eval(exampleName); GlobalConstants.setDummyMode(false);
fname = fullfile(outDir,'ExampleClosedModel.java');
fid = fopen(fname,'a');
fprintf(fid,'\n\t@Test\n');
fprintf(fid,'\tpublic void test_%s() throws IllegalAccessException {\n',exampleName);
QN2JAVA(model, exampleName, fid, false);
printExample(fid, exampleName, model, solver, AvgTable);
fprintf(fid, '\t}\n');
fclose(fid);
%%
GlobalConstants.setDummyMode(true); exampleName='cqn_threeclass_hyperl'; eval(exampleName); GlobalConstants.setDummyMode(false);
fname = fullfile(outDir,'ExampleClosedModel.java');
fid = fopen(fname,'a');
fprintf(fid,'\n\t@Test\n');
fprintf(fid,'\tpublic void test_%s() throws IllegalAccessException {\n',exampleName);
QN2JAVA(model, exampleName, fid, false);
printExample(fid, exampleName, model, solver, AvgTable);
fprintf(fid, '\t}\n');
fclose(fid);
%%
fname = fullfile(outDir,'ExampleClosedModel.java');
fid = fopen(fname,'a');
fprintf(fid, '}\n');
if fid~=1
    fclose(fid);
end

function printExample(fid, exampleName, model, solver, avgTable)
%for s=1:length(solver)
for s = 5
    solver{s} = SolverMVA(model,'qdlin','verbose',false);
    avgTable{s} = solver{s}.getAvgTable;
    fprintf(fid,'\t\t%s solver = new %s(model);\n\n',solver{s}.name,solver{s}.name);
    fprintf(fid,'\t\tList<List<Double>> avgTable = solver.getAvgTable();\n\n');
    fprintf(fid,'\t\tList<Double> QLen = avgTable.get(0);\n');

    QLen = avgTable{s}.QLen;
    for i=1:length(QLen)
        fprintf(fid,'\t\tassertEquals(%.16f, QLen.get(%d));\n',QLen(i),i-1);
    end
    fprintf(fid, '\n');

    fprintf(fid,'\t\tList<Double> Util = avgTable.get(1);\n');
    Util = avgTable{s}.Util;
    for i=1:length(Util)
        fprintf(fid,'\t\tassertEquals(%.16f, Util.get(%d));\n',Util(i),i-1);
    end
    fprintf(fid, '\n');

    fprintf(fid,'\t\tList<Double> RespT = avgTable.get(2);\n');
    RespT = avgTable{s}.RespT;
    for i=1:length(RespT)
        fprintf(fid,'\t\tassertEquals(%.16f, RespT.get(%d));\n',RespT(i),i-1);
    end
    fprintf(fid, '\n');

    fprintf(fid,'\t\tList<Double> ResidT = avgTable.get(3);\n');
    ResidT = avgTable{s}.ResidT;
    for i=1:length(ResidT)
        fprintf(fid,'\t\tassertEquals(%.16f, ResidT.get(%d));\n',ResidT(i),i-1);
    end
    fprintf(fid, '\n');

    %     fprintf(fid,'\t\tList<Double> ArvR = avgTable.get(4);\n');
    %     ArvR = AvgTable{s}.ArvR;
    %     for i=1:length(ArvR)
    %         fprintf(fid,'\t\tassertEquals(%.16f, ArvR.get(%d));\n',ArvR(i),i-1);
    %     end
    %     fprintf(fid, '\n');

    fprintf(fid,'\t\tList<Double> Tput = avgTable.get(4);\n');
    %fprintf(fid,'\t\tList<Double> Tput = avgTable.get(5);\n');
    Tput = avgTable{s}.Tput;
    for i=1:length(Tput)
        fprintf(fid,'\t\tassertEquals(%.16f, Tput.get(%d));\n',Tput(i),i-1);
    end
end
end