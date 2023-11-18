import org.apache.commons.cli.*;
import java.util.*;
import java.io.*;
import java.net.*;
import org.java_websocket.client.WebSocketClient;
import org.java_websocket.drafts.Draft;
import org.java_websocket.handshake.ServerHandshake;

/**
 * This example demonstrates how to create a websocket connection to a server. Only the most
 * important callbacks are overloaded.
 */
public class lineclient extends WebSocketClient {
  
 public String message;

  public lineclient(URI serverURI, String msg) {
    super(serverURI);
    this.message = msg;
  }

  @Override
  public void onOpen(ServerHandshake handshakedata) { 
  	send(this.message);
  }

  @Override
  public void onMessage(String message) { System.out.println(message); }

  @Override
  public void onClose(int code, String reason, boolean remote) { }

  @Override
  public void onError(Exception ex) { ex.printStackTrace(); }

  public static void main(String[] args) throws URISyntaxException, ParseException {
        Options options = new Options();
	options.addOption("h", false, "Help");
	options.addOption("m", false, "Max number of requests");
	options.addOption("s", false, "Solver");
	options.addOption("a", false, "Analysis");
	options.addOption("f", false, "File");
	options.addOption("v", false, "Verbosity");
	options.addOption("i", false, "Input");
	options.addOption("o", false, "Output");
	options.addOption("d", false, "Seed");

	CommandLineParser parser = new DefaultParser();
	CommandLine cmd = parser.parse(options, args);

	if (cmd.hasOption("h")) {
		System.out.println("LINE Solver - Client Interface");
		System.out.println("Copyright (c) 2012-2021, QORE group, Imperial College London");
		System.out.println("-----------------------------------------------------------------------");
		System.out.println("Usage: cat MODEL | java -jar linecli.jar SERVER_IP SERVER_PORT\n");
		System.out.println("Example: cat mymodel.jsimg | java -jar linecli.jar 192.168.0.1 5863 ");
	}
	else {
		
		String message = String.join(",", args);
		Scanner sc = new Scanner(System.in);
        	for(int i = 1; sc.hasNext()== true; i++){
            		message = message + "\n" + sc.nextLine();
        	}
    		lineclient c = new lineclient(new URI("ws://"+args[0]+":"+args[1]), message); 
    		c.connect();
  	}
 }
}
