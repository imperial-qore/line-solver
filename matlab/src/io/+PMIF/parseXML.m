function sn = parseXML(filename, verbose)
% QN = PARSEXML(FILENAME, VERBOSE)

% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.DocumentBuilder;
import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import org.w3c.dom.Node;
import org.w3c.dom.Element;
import java.io.File;


dbFactory = DocumentBuilderFactory.newInstance();
dBuilder = dbFactory.newDocumentBuilder();
try
    doc = dBuilder.parse(filename);
catch exception %java.io.FileNotFoundException
    if ~exist(filename, 'file')
        line_printf(['Error: Input file ', filename,' not found']);
        sn = [];
        return;
    else
        rethrow(exception);
    end
end

doc.getDocumentElement().normalize();
rootElem = char(doc.getDocumentElement().getNodeName());
if verbose > 0
    line_printf(['Root element:', rootElem ] );
end

if ~strcmp(rootElem, 'QueueingNetworkModel')
    line_printf(['Error: Input file ', filename,' does not have a proper PMIF model']);
    line_printf(['Root element:', rootElem ] );
    sn = [];
    return;
end


% read servers
serverList = doc.getElementsByTagName('Server');
servers = [];
for i = 0:serverList.getLength()-1
    serverNode = serverList.item(i);
    if serverNode.getNodeType() == Node.ELEMENT_NODE
        serverElement = serverNode;
        name = char(serverElement.getAttribute('Name'));
        quantity = str2num(char(serverElement.getAttribute('Quantity')));
        scheduling = char(serverElement.getAttribute('SchedulingPolicy'));
        tempServer = PMIF.server(name, quantity, scheduling);
        servers = [servers; tempServer];
    end
end

workUnitServerList = doc.getElementsByTagName('WorkUnitServer');
workUnitServers = [];
for i = 0:workUnitServerList.getLength()-1
    serverNode = workUnitServerList.item(i);
    if serverNode.getNodeType() == Node.ELEMENT_NODE
        serverElement = serverNode;
        name = char(serverElement.getAttribute('Name'));
        quantity = str2num(char(serverElement.getAttribute('Quantity')));
        scheduling = char(serverElement.getAttribute('SchedulingPolicy'));
        serviceTime = str2num(char(serverElement.getAttribute('ServiceTime')));
        timeUnits = char(serverElement.getAttribute('TimeUnits'));
        tempServer = PMIF.workUnitServer(name, quantity, scheduling,serviceTime,timeUnits);
        workUnitServers = [workUnitServers; tempServer];
    end
end

closedWorkloadList = doc.getElementsByTagName('ClosedWorkload');
closedWorkloads = [];
for i = 0:closedWorkloadList.getLength()-1
    workloadNode = closedWorkloadList.item(i);
    if workloadNode.getNodeType() == Node.ELEMENT_NODE
        workloadElement = workloadNode;
        name = char(workloadElement.getAttribute('WorkloadName'));
        numberJobs = str2num(char(workloadElement.getAttribute('NumberOfJobs')));
        thinkTime = str2num(char(workloadElement.getAttribute('ThinkTime')));
        thinkDevice = char(workloadElement.getAttribute('ThinkDevice'));
        timeUnits = char(workloadElement.getAttribute('TimeUnits'));
        tempWorkload = PMIF.closedWorkload(name, numberJobs, thinkTime, thinkDevice, timeUnits);
        tempWorkload = addTransits(tempWorkload, workloadNode);
        closedWorkloads = [closedWorkloads; tempWorkload];
    end
end

demandServiceRequestList = doc.getElementsByTagName('DemandServiceRequest');
demandServiceRequests = [];
for i = 0:demandServiceRequestList.getLength()-1
    demandNode = demandServiceRequestList.item(i);
    if demandNode.getNodeType() == Node.ELEMENT_NODE
        demandElement = demandNode;
        workloadName = char(demandElement.getAttribute('WorkloadName'));
        serverID = char(demandElement.getAttribute('ServerID'));
        serviceDemand = str2num(char(demandElement.getAttribute('ServiceDemand')));
        numberVisits = str2num(char(demandElement.getAttribute('NumberOfVisits')));
        timeUnits = char(workloadElement.getAttribute('TimeUnits'));
        tempDemand = PMIF.demandServiceRequest(workloadName, serverID, serviceDemand, numberVisits, timeUnits);
        tempDemand = addTransits(tempDemand, demandNode);
        demandServiceRequests = [demandServiceRequests; tempDemand];
    end
end

workUnitServiceRequestList = doc.getElementsByTagName('WorkUnitServiceRequest');
workUnitServiceRequests = [];
for i = 0:workUnitServiceRequestList.getLength()-1
    demandNode = workUnitServiceRequestList.item(i);
    if demandNode.getNodeType() == Node.ELEMENT_NODE
        demandElement = demandNode;
        workloadName = char(demandElement.getAttribute('WorkloadName'));
        serverID = char(demandElement.getAttribute('ServerID'));
        numberVisits = str2num(char(demandElement.getAttribute('NumberOfVisits')));
        tempDemand = PMIF.workUnitServiceRequest(workloadName, serverID, numberVisits);
        tempDemand = addTransits(tempDemand, demandNode);
        workUnitServiceRequests = [workUnitServiceRequests; tempDemand];
    end
end


timeServiceRequestList = doc.getElementsByTagName('TimeServiceRequest');
timeServiceRequests = [];
for i = 0:timeServiceRequestList.getLength()-1
    demandNode = timeServiceRequestList.item(i);
    if demandNode.getNodeType() == Node.ELEMENT_NODE
        demandElement = demandNode;
        workloadName = char(demandElement.getAttribute('WorkloadName'));
        serverID = char(demandElement.getAttribute('ServerID'));
        serviceDemand = str2num(char(demandElement.getAttribute('ServiceDemand')));
        timeUnits = char(workloadElement.getAttribute('TimeUnits'));
        tempDemand = timeServiceRequest(workloadName, serverID, serviceDemand, timeUnits);
        tempDemand = addTransits(tempDemand, demandNode);
        timeServiceRequests = [timeServiceRequests; tempDemand];
    end
end


sn = PMIF.PMIF(servers, workUnitServers, closedWorkloads, ...
    demandServiceRequests, workUnitServiceRequests, timeServiceRequests);


end

function tempObj = addTransits(tempObj, node)
% TEMPOBJ = ADDTRANSITS(TEMPOBJ, NODE)

import org.w3c.dom.Node;

% get transits
transitList = node.getElementsByTagName('Transit');
for j = 0:transitList.getLength()-1
    transitNode = transitList.item(j);
    if transitNode.getNodeType() == Node.ELEMENT_NODE
        transitElement = transitNode;
        dest = char(transitElement.getAttribute('To'));
        prob = str2num(char(transitElement.getAttribute('Probability')));
        tempObj = tempObj.addTransit(dest, prob);
    end
end
end
