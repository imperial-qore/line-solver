package jline.lang.constant;

import java.io.Serializable;

/**
 *  Constants for specifying a join strategy
 */
public enum JoinStrategy implements Serializable {
	STD,
	PARTIAL,
	Quorum,
	Guard
}
