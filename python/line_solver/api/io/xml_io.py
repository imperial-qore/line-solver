"""
XML I/O for LINE Network models.

This module provides XML import/export functionality for Network models,
enabling model serialization, sharing, and interoperability.

Port from:
    - matlab/src/io/NetworkXMLIO.m
"""

import xml.etree.ElementTree as ET
from xml.dom import minidom
from typing import Optional, Dict, Any, List, Tuple
from dataclasses import dataclass
from pathlib import Path
import numpy as np


@dataclass
class NetworkXMLInfo:
    """Information about a Network XML file."""
    name: str
    version: str
    node_count: int
    class_count: int
    connection_count: int
    node_types: Dict[str, int]


class NetworkXMLIO:
    """
    XML import/export for LINE Network models.

    Provides functionality to save and load Network models to/from XML format.
    Supports class switching through ClassSwitch nodes and cross-class routing.

    Usage:
        # Export a network to XML
        NetworkXMLIO.export_to_xml(model, 'mymodel.xml')

        # Import a network from XML
        model = NetworkXMLIO.import_from_xml('mymodel.xml')

        # Get XML file info without full import
        info = NetworkXMLIO.get_xml_info('mymodel.xml')

    References:
        MATLAB: matlab/src/io/NetworkXMLIO.m
    """

    XML_VERSION = "1.0"
    NAMESPACE = "http://line.imperial.ac.uk/network"

    @staticmethod
    def export_to_xml(network: Any, filename: str) -> None:
        """
        Export a Network model to XML file.

        Args:
            network: Network model to export
            filename: Output XML file path

        Raises:
            ValueError: If network is not a valid Network object
            IOError: If export fails
        """
        try:
            # Get network structure
            if hasattr(network, 'getStruct'):
                sn = network.getStruct()
            else:
                sn = network

            # Create root element
            root = ET.Element('network')
            root.set('name', getattr(sn, 'name', 'model'))
            root.set('version', NetworkXMLIO.XML_VERSION)
            root.set('xmlns', NetworkXMLIO.NAMESPACE)

            # Add metadata
            metadata = ET.SubElement(root, 'metadata')
            ET.SubElement(metadata, 'nstations').text = str(sn.nstations)
            ET.SubElement(metadata, 'nclasses').text = str(sn.nclasses)
            ET.SubElement(metadata, 'nnodes').text = str(sn.nnodes)

            # Add nodes
            nodes_elem = ET.SubElement(root, 'nodes')
            for i in range(sn.nnodes):
                node_elem = ET.SubElement(nodes_elem, 'node')
                node_elem.set('id', str(i))
                node_elem.set('name', sn.nodenames[i] if hasattr(sn, 'nodenames') else f'Node{i}')

                # Node type
                if hasattr(sn, 'nodetype'):
                    node_type = sn.nodetype[i]
                    node_elem.set('type', str(node_type.name if hasattr(node_type, 'name') else node_type))

                # Scheduling strategy (for stations)
                if hasattr(sn, 'nodeToStation') and hasattr(sn, 'sched'):
                    try:
                        ist = sn.nodeToStation[i]
                        if ist is not None and ist < len(sn.sched):
                            sched = sn.sched[ist]
                            node_elem.set('scheduling', str(sched.name if hasattr(sched, 'name') else sched))
                    except (IndexError, TypeError):
                        pass

                # Number of servers
                if hasattr(sn, 'nservers'):
                    try:
                        ist = sn.nodeToStation[i] if hasattr(sn, 'nodeToStation') else i
                        if ist is not None and ist < len(sn.nservers):
                            node_elem.set('servers', str(int(sn.nservers[ist])))
                    except (IndexError, TypeError):
                        pass

            # Add classes
            classes_elem = ET.SubElement(root, 'classes')
            for k in range(sn.nclasses):
                class_elem = ET.SubElement(classes_elem, 'class')
                class_elem.set('id', str(k))
                class_elem.set('name', sn.classnames[k] if hasattr(sn, 'classnames') else f'Class{k}')

                # Population
                if hasattr(sn, 'njobs'):
                    njobs = sn.njobs[k]
                    class_elem.set('population', 'inf' if np.isinf(njobs) else str(int(njobs)))
                    class_elem.set('type', 'open' if np.isinf(njobs) else 'closed')

                # Priority
                if hasattr(sn, 'classprio'):
                    class_elem.set('priority', str(sn.classprio[k]))

                # Reference station
                if hasattr(sn, 'refstat') and not np.isinf(sn.njobs[k]):
                    class_elem.set('refstation', str(int(sn.refstat[k])))

            # Add service/arrival processes
            processes_elem = ET.SubElement(root, 'processes')
            if hasattr(sn, 'proc') and sn.proc is not None:
                for ist in range(sn.nstations):
                    for k in range(sn.nclasses):
                        if ist < len(sn.proc) and sn.proc[ist] is not None:
                            if k < len(sn.proc[ist]) and sn.proc[ist][k] is not None:
                                proc_elem = ET.SubElement(processes_elem, 'process')
                                proc_elem.set('station', str(ist))
                                proc_elem.set('class', str(k))

                                # Store rate
                                if hasattr(sn, 'rates') and sn.rates is not None:
                                    rate = sn.rates[ist, k]
                                    if not np.isnan(rate):
                                        proc_elem.set('rate', str(rate))

                                # Store SCV
                                if hasattr(sn, 'scv') and sn.scv is not None:
                                    scv = sn.scv[ist, k]
                                    if not np.isnan(scv):
                                        proc_elem.set('scv', str(scv))

            # Add routing (connections)
            routing_elem = ET.SubElement(root, 'routing')
            if hasattr(sn, 'rtnodes') and sn.rtnodes is not None:
                rtnodes = sn.rtnodes
                for i in range(sn.nnodes):
                    for k in range(sn.nclasses):
                        for j in range(sn.nnodes):
                            for c in range(sn.nclasses):
                                idx_from = i * sn.nclasses + k
                                idx_to = j * sn.nclasses + c
                                if idx_from < rtnodes.shape[0] and idx_to < rtnodes.shape[1]:
                                    prob = rtnodes[idx_from, idx_to]
                                    if prob > 0:
                                        link_elem = ET.SubElement(routing_elem, 'link')
                                        link_elem.set('from_node', str(i))
                                        link_elem.set('from_class', str(k))
                                        link_elem.set('to_node', str(j))
                                        link_elem.set('to_class', str(c))
                                        link_elem.set('probability', str(prob))

            # Write to file with pretty formatting
            xml_str = ET.tostring(root, encoding='unicode')
            dom = minidom.parseString(xml_str)
            pretty_xml = dom.toprettyxml(indent='  ')

            with open(filename, 'w', encoding='utf-8') as f:
                f.write(pretty_xml)

            print(f"Network successfully exported to: {filename}")

        except Exception as e:
            raise IOError(f"Failed to export network: {e}")

    @staticmethod
    def import_from_xml(filename: str) -> Any:
        """
        Import a Network model from XML file.

        Args:
            filename: Input XML file path

        Returns:
            Network model (as NetworkStruct-like object)

        Raises:
            IOError: If import fails
            ValueError: If XML structure is invalid
        """
        try:
            tree = ET.parse(filename)
            root = tree.getroot()

            # Validate root element
            if root.tag != 'network':
                raise ValueError("Root element must be <network>")

            # Create network structure dictionary
            network_data = {
                'name': root.get('name', 'imported_model'),
                'version': root.get('version', '1.0'),
            }

            # Parse metadata
            metadata = root.find('metadata')
            if metadata is not None:
                network_data['nstations'] = int(metadata.findtext('nstations', '0'))
                network_data['nclasses'] = int(metadata.findtext('nclasses', '0'))
                network_data['nnodes'] = int(metadata.findtext('nnodes', '0'))

            # Parse nodes
            nodes_elem = root.find('nodes')
            nodes = []
            if nodes_elem is not None:
                for node_elem in nodes_elem.findall('node'):
                    node_data = {
                        'id': int(node_elem.get('id', '0')),
                        'name': node_elem.get('name', ''),
                        'type': node_elem.get('type', 'Queue'),
                        'scheduling': node_elem.get('scheduling', 'FCFS'),
                        'servers': int(node_elem.get('servers', '1')),
                    }
                    nodes.append(node_data)
            network_data['nodes'] = nodes

            # Parse classes
            classes_elem = root.find('classes')
            classes = []
            if classes_elem is not None:
                for class_elem in classes_elem.findall('class'):
                    pop_str = class_elem.get('population', '0')
                    population = np.inf if pop_str == 'inf' else int(pop_str)
                    class_data = {
                        'id': int(class_elem.get('id', '0')),
                        'name': class_elem.get('name', ''),
                        'population': population,
                        'type': class_elem.get('type', 'closed'),
                        'priority': int(class_elem.get('priority', '0')),
                        'refstation': int(class_elem.get('refstation', '0')),
                    }
                    classes.append(class_data)
            network_data['classes'] = classes

            # Parse processes
            processes_elem = root.find('processes')
            processes = []
            if processes_elem is not None:
                for proc_elem in processes_elem.findall('process'):
                    proc_data = {
                        'station': int(proc_elem.get('station', '0')),
                        'class': int(proc_elem.get('class', '0')),
                        'rate': float(proc_elem.get('rate', '1.0')),
                        'scv': float(proc_elem.get('scv', '1.0')),
                    }
                    processes.append(proc_data)
            network_data['processes'] = processes

            # Parse routing
            routing_elem = root.find('routing')
            links = []
            if routing_elem is not None:
                for link_elem in routing_elem.findall('link'):
                    link_data = {
                        'from_node': int(link_elem.get('from_node', '0')),
                        'from_class': int(link_elem.get('from_class', '0')),
                        'to_node': int(link_elem.get('to_node', '0')),
                        'to_class': int(link_elem.get('to_class', '0')),
                        'probability': float(link_elem.get('probability', '0.0')),
                    }
                    links.append(link_data)
            network_data['routing'] = links

            print(f"Network successfully imported from: {filename}")
            return network_data

        except ET.ParseError as e:
            raise IOError(f"XML parsing error: {e}")
        except Exception as e:
            raise IOError(f"Failed to import network: {e}")

    @staticmethod
    def validate_xml(filename: str) -> bool:
        """
        Validate XML file structure without full import.

        Args:
            filename: XML file path to validate

        Returns:
            True if valid

        Raises:
            ValueError: If validation fails
        """
        try:
            tree = ET.parse(filename)
            root = tree.getroot()

            if root.tag != 'network':
                raise ValueError("Root element must be <network>")

            # Check required elements
            required = ['nodes', 'classes']
            for elem_name in required:
                if root.find(elem_name) is None:
                    raise ValueError(f"Missing required element: <{elem_name}>")

            print(f"XML file is valid: {filename}")
            return True

        except ET.ParseError as e:
            raise ValueError(f"XML validation failed: {e}")

    @staticmethod
    def get_xml_info(filename: str) -> NetworkXMLInfo:
        """
        Get basic information about XML file without full import.

        Args:
            filename: XML file path

        Returns:
            NetworkXMLInfo with network information
        """
        try:
            tree = ET.parse(filename)
            root = tree.getroot()

            name = root.get('name', '')
            version = root.get('version', '')

            # Count nodes
            nodes_elem = root.find('nodes')
            node_count = len(nodes_elem.findall('node')) if nodes_elem is not None else 0

            # Count classes
            classes_elem = root.find('classes')
            class_count = len(classes_elem.findall('class')) if classes_elem is not None else 0

            # Count connections
            routing_elem = root.find('routing')
            connection_count = len(routing_elem.findall('link')) if routing_elem is not None else 0

            # Get node types
            node_types: Dict[str, int] = {}
            if nodes_elem is not None:
                for node_elem in nodes_elem.findall('node'):
                    node_type = node_elem.get('type', 'Unknown')
                    node_types[node_type] = node_types.get(node_type, 0) + 1

            return NetworkXMLInfo(
                name=name,
                version=version,
                node_count=node_count,
                class_count=class_count,
                connection_count=connection_count,
                node_types=node_types
            )

        except Exception as e:
            raise IOError(f"Failed to get XML info: {e}")


# Convenience functions
def export_to_xml(network: Any, filename: str) -> None:
    """Export a Network model to XML file."""
    NetworkXMLIO.export_to_xml(network, filename)


def import_from_xml(filename: str) -> Any:
    """Import a Network model from XML file."""
    return NetworkXMLIO.import_from_xml(filename)


def validate_xml(filename: str) -> bool:
    """Validate XML file structure."""
    return NetworkXMLIO.validate_xml(filename)


def get_xml_info(filename: str) -> NetworkXMLInfo:
    """Get basic information about XML file."""
    return NetworkXMLIO.get_xml_info(filename)


__all__ = [
    'NetworkXMLIO',
    'NetworkXMLInfo',
    'export_to_xml',
    'import_from_xml',
    'validate_xml',
    'get_xml_info',
]
