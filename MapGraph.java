package roadgraph;


import java.util.*;
import java.util.function.Consumer;

import geography.GeographicPoint;
import sun.awt.image.ImageWatched;
import util.GraphLoader;

/**
 * @author Noah Garner
 *
 * A class which represents a graph of geographic locations
 * Nodes in the graph are intersections between
 */
public class MapGraph {
    // Maintain both nodes and edges as you will need to
    // be able to look up nodes by lat/lon or by roads
    // that contain those nodes.
    private HashMap<GeographicPoint,MapNode> pointNodeMap; //geo point takes a long/lat as key, and the value is node
    private HashSet<MapEdge> edges;

    /**
     * Create a new empty MapGraph
     */

    public MapGraph()
    {
        pointNodeMap = new HashMap<>();
        edges = new HashSet<>();
    }

    /**
     * Get the number of vertices (road intersections) in the graph
     * @return The number of vertices in the graph.
     */
    public int getNumVertices()
    {
        return pointNodeMap.values().size();
    }

    /**
     * Return the intersections, which are the vertices in this graph.
     * @return The vertices in this graph as GeographicPoints
     */
    public Set<GeographicPoint> getVertices()
    {
        return pointNodeMap.keySet();
    }

    /**
     * Get the number of road segments in the graph
     * @return The number of edges in the graph.
     */
    public int getNumEdges()
    {
        return edges.size();
    }



    /** Add a node corresponding to an intersection at a Geographic Point
     * If the location is already in the graph or null, this method does
     * not change the graph.
     * @param location  The location of the intersection
     * @return true if a node was added, false if it was not (the node
     * was already in the graph, or the parameter is null).
     */
    public boolean addVertex(GeographicPoint location)
    {
        if (location == null) {
            return false;
        }
        MapNode n = pointNodeMap.get(location);
        if (n == null) {
            n = new MapNode(location);
            pointNodeMap.put(location, n);
            return true;
        }
        else {
            return false;
        }
    }

    /**
     * Adds a directed edge to the graph from pt1 to pt2.
     * Precondition: Both GeographicPoints have already been added to the graph
     * @param from The starting point of the edge
     * @param to The ending point of the edge
     * @param roadName The name of the road
     * @param roadType The type of the road
     * @param length The length of the road, in km
     * @throws IllegalArgumentException If the points have not already been
     *   added as nodes to the graph, if any of the arguments is null,
     *   or if the length is less than 0.
     */

    public void addEdge(GeographicPoint from, GeographicPoint to, String roadName,
                        String roadType, double length) throws IllegalArgumentException {

        MapNode n1 = pointNodeMap.get(from);
        MapNode n2 = pointNodeMap.get(to);

        // check nodes are valid
        if (n1 == null)
            throw new NullPointerException("addEdge: pt1:"+from+"is not in graph");
        if (n2 == null)
            throw new NullPointerException("addEdge: pt2:"+to+"is not in graph");

        MapEdge edge = new MapEdge(roadName, roadType, n1, n2, length);
        edges.add(edge);
        n1.addEdge(edge);

    }

    /**
     * Get a set of neighbor nodes from a mapNode
     * @param node  The node to get the neighbors from
     * @return A set containing the MapNode objects that are the neighbors
     * 	of node
     */
    private Set<MapNode> getNeighbors(MapNode node) {
        return node.getNeighbors();
    }

    /** Find the path from start to goal using breadth first search
     *
     * @param start The starting location
     * @param goal The goal location
     * @return The list of intersections that form the shortest (unweighted)
     *   path from start to goal (including both start and goal).
     */
    public List<GeographicPoint> bfs(GeographicPoint start, GeographicPoint goal) {
        // Dummy variable for calling the search algorithms
        Consumer<GeographicPoint> temp = (x) -> {};
        return bfs(start, goal, temp);
    }

    public List<GeographicPoint> bfs(GeographicPoint start,
                                     GeographicPoint goal,
                                     Consumer<GeographicPoint> nodeSearched) {
        if (!checkPoints(start, goal) ) return null;
        //initialize everything we will use, we begin by creating start/goal nodes
        MapNode startNode = pointNodeMap.get(start);
        MapNode endNode = pointNodeMap.get(goal);

        // setup to begin BFS
        HashMap<MapNode,MapNode> parentMap = new HashMap<>();
        Queue<MapNode> toExplore = new LinkedList<>();
        HashSet<MapNode> visited = new HashSet<>();

        toExplore.add(startNode);
        MapNode currNode = null;

        while (!toExplore.isEmpty()) {
            currNode = toExplore.remove();

            // hook for visualization
            nodeSearched.accept(currNode.getLocation());

            if (currNode.equals(endNode)) break; //we found goal!
            Set<MapNode> neighbors = getNeighbors(currNode);
            for (MapNode neighbor : neighbors) {
                if (!visited.contains(neighbor)) {
                    visited.add(neighbor);          //add to visited list
                    parentMap.put(neighbor, currNode);  //the value is the currNode, the neighbor is the key
                                                    //because guess what, its a parentMap, so we should have the value of the key be a parent of the key.
                                                    //in the context of the map, since neighbor is being looked at, neighbor's value in parentMap must be
                                                    //the previous node from hence it came. This concept is whats been taking me so fucking long to understand
                                                    //but I think now I finally understand it
                    toExplore.add(neighbor);        //queue up the neighbors of the currNode!
                }
            }
        }
        if (!currNode.equals(endNode)) {
            System.out.println("No path found from " +start+ " to " + goal);
            return null;
        }
        // Reconstruct the parent path
        return  reconstructPath(parentMap, startNode, endNode);

    }


    /** Find the path from start to goal using Dijkstra's algorithm
     *
     * @param start The starting location
     * @param goal The goal location
     * @return The list of intersections that form the shortest path from
     *   start to goal (including both start and goal).
     */
    public List<GeographicPoint> dijkstra(GeographicPoint start, GeographicPoint goal) {
        // Dummy variable for calling the search algorithms
        // You do not need to change this method.
        Consumer<GeographicPoint> temp = (x) -> {};
        return dijkstra(start, goal, temp);
    }


    public List<GeographicPoint> dijkstra(GeographicPoint start,
                                          GeographicPoint goal, Consumer<GeographicPoint> nodeSearched) {

        if (!checkPoints(start, goal)) return null;
        //initialize everything we will use, we begin by creating start/goal nodes
        MapNode startNode = pointNodeMap.get(start);
        MapNode endNode = pointNodeMap.get(goal);

        HashSet<MapNode> visited = new HashSet<>(); //track already visited
        HashMap<MapNode, MapNode> parentMap = new HashMap<>(); //track each nodes parent
        Queue<MapNode> pq = new PriorityQueue<>(); //this queue implements comparator

        //so we can re-initialize if needed
        for (MapNode n : pointNodeMap.values()){
            n.setDistance(Double.POSITIVE_INFINITY);
            n.setActualDistance(Double.POSITIVE_INFINITY);
        }

        startNode.setActualDistance(0.0); //initialize our start with 0 weight
        startNode.setDistance(0.0);
        pq.add(startNode); // we want to add our start to the queue, this is when compareTo() is called
        MapNode currNode = null;
        while(!pq.isEmpty()) {
            //because MapNode implements comparator, and this is a priority queue,
            //we know that it removes the node with the smallest distance
            currNode = pq.remove(); //take off queue, this takes the top element off, and puts the bottom element in the heap at top, and bubbles it down
            nodeSearched.accept(currNode.getLocation()); //for visualization

            if(!visited.contains(currNode)) {
                visited.add(currNode);

                if(currNode.equals(endNode))
                    break;
                Set<MapEdge> currNodeEdges = currNode.getEdges(); //

                for (MapEdge edge : currNodeEdges) { //for each neighbor in neighbors of currNode
                    MapNode aNeighbor = edge.getOtherNode(currNode); //use this edge, and currNode to get the neighbor at other end of the edge
                    //if this neighbor has not been visited
                    if (!visited.contains(aNeighbor)) {
                        //update aNeighbor's distance:
                        //with the currNode's weighted distance, and the edge's length to the neighbor being looked at
                        double newDistance = currNode.getDistance() + edge.getLength();
                        //if this distance is less than the neighbor's, we update iaNeighbor's distance, and add it to parent map
                        //so we can reach it, and add aNeighbor to the queue
                        if (newDistance < aNeighbor.getDistance()) {
                            aNeighbor.setDistance(newDistance);
                            aNeighbor.setActualDistance(newDistance);
                            parentMap.put(aNeighbor, currNode);
                            pq.add(aNeighbor); //arbitrary added to queue data structure
                        }
                    }
                }
            }
        }
        if (!currNode.equals(endNode)) {
            System.out.println("No path found from " +start+ " to " + goal);
            return null;
        }
        System.out.println("Nodes Visited Dijkstra: "+visited.size());
        //return the reconstructed path
        return reconstructPath(parentMap, startNode, endNode);
    }

    /** Find the path from start to goal using A-Star search
     *
     * @param start The starting location
     * @param goal The goal location
     * @return The list of intersections that form the shortest path from
     *   start to goal (including both start and goal).
     */
    public List<GeographicPoint> aStarSearch(GeographicPoint start, GeographicPoint goal) {
        Consumer<GeographicPoint> temp = (x) -> {
            // Dummy variable for calling the search algorithms
        };
        return aStarSearch(start, goal, temp);
    }


    public List<GeographicPoint> aStarSearch(GeographicPoint start,
                                             GeographicPoint goal, Consumer<GeographicPoint> nodeSearched) {
        if (!checkPoints(start, goal)) return null; //send to our check function
        //check went through, lets make our MapNodes out of our geo points so we can traverse from the s -> g
        MapNode startNode = pointNodeMap.get(start);
        MapNode endNode =pointNodeMap.get(goal);
        //initialize data structures'
        PriorityQueue<MapNode> pq = new PriorityQueue<>(); //implemented as a min-heap, the compareTo method is called with each addition
        HashMap<MapNode,MapNode> parentMap = new HashMap<>();
        HashSet<MapNode> visited = new HashSet<>();
        //so we can re-initialize if needed, all nodes in a weighted graph using Dikstra and A* must start at infinity
        for (MapNode nodeToInitialize : pointNodeMap.values()){
            nodeToInitialize.setDistance(Double.POSITIVE_INFINITY);
            nodeToInitialize.setActualDistance(Double.POSITIVE_INFINITY);
        }

        //initialize start of our procedure
        startNode.setActualDistance(0.0);//remember the start has no weight!
        startNode.setDistance(0.0);
        pq.add(startNode); //we want to add the start  to queue to begin algorithm
        MapNode currNode = null; //establish a node that will be used for traversing the graph
        //while there are nodesi n queue
        while(!pq.isEmpty()) {
            currNode = pq.remove();
            nodeSearched.accept(currNode.getLocation()); //for visualization
            //if it has not been visited
            if(!visited.contains(currNode)) {
                visited.add(currNode);
                //if it is the end, break out and reconstruct the path, otherwise continue to creating set
                if (currNode.equals(endNode)) //to stop the while loop
                    break; //done, leave function
                visitAStar(visited, currNode, goal, parentMap, pq);
            }
        }
        if (!currNode.equals(endNode)) {
            System.out.println("No path found from " +start+ " to " + goal);
            return null;
        }
        System.out.println("Nodes Visited A Star: "+visited.size());
        //return the reconstructed path
        return reconstructPath(parentMap, startNode, endNode);
    }
    //refactored out
    private List<GeographicPoint> reconstructPath(HashMap<MapNode,
            MapNode> parentMap, MapNode start, MapNode goal) {
        LinkedList<GeographicPoint> path = new LinkedList<GeographicPoint>();
        MapNode current = goal;

        while (!current.equals(start)) {
            path.addFirst(current.getLocation());
            current = parentMap.get(current);
        }
        // add start
        path.addFirst(start.getLocation());
        return path;
    }

    public double aStarSearchExtension(GeographicPoint start, GeographicPoint goal) {
        // Dummy variable for calling the search algorithms
        Consumer<GeographicPoint> temp = (x) -> {};
        return aStarSearchExtension(start, goal, temp);
    }

    public double aStarSearchExtension(GeographicPoint start,
                                             GeographicPoint goal, Consumer<GeographicPoint> nodeSearched) {
        if (!checkPoints(start, goal)) return 0; //send to our check function
        //check went through, lets make our MapNodes out of our geo points so we can traverse from the s -> g
        MapNode startNode = pointNodeMap.get(start);
        MapNode endNode =pointNodeMap.get(goal);
        //initialize data structures'
        PriorityQueue<MapNode> pq = new PriorityQueue<>(); //implemented as a min-heap, the compareTo method is called with each addition
        HashMap<MapNode,MapNode> parentMap = new HashMap<>();
        HashSet<MapNode> visited = new HashSet<>();
        //so we can re-initialize if needed, all nodes in a weighted graph using Dikstra and A* must start at infinity
        //initialize start of our procedure
        for (MapNode nodeToInitialize : pointNodeMap.values()){
            nodeToInitialize.setDistance(Double.POSITIVE_INFINITY);
            nodeToInitialize.setActualDistance(Double.POSITIVE_INFINITY);
        }
        startNode.setActualDistance(0.0);//remember the start has no weight!
        startNode.setDistance(0.0);
        pq.add(startNode); //we want to ad
        // d the start  to queue to begin algorithm
        MapNode currNode = null; //establish a node that will be used for traversing the graph
        //while there are nodesi n queue
        while(!pq.isEmpty()) {
            currNode = pq.remove();
            nodeSearched.accept(currNode.getLocation()); //for visualization
            //if it has not been visited
            if(!visited.contains(currNode)) {
                visited.add(currNode);
                //if it is the end, break out and reconstruct the path, otherwise continue to creating set
                if (currNode.equals(endNode)) //to stop the while loop
                    break; //done, leave function
                visitAStar(visited, currNode, goal, parentMap, pq);
            }
        }
        if (!currNode.equals(endNode)) {
            return 0;
        }
        //return the reconstructed path
        return reconstructPathExtension(parentMap, startNode, endNode);
    }
    //helper
    //Refactoring a big piece of a function into this smaller function:
    // Steps
    //helper methods for both of my A* implementations
    // 1. copy paste code, give it a function name of meaning
    // 2. Enter function parameters as needed to satisfy errors
    // 3. Call the new function in the old function. Done
    //refactored out function for dijkstra algorithm
    public void visitAStar(HashSet<MapNode> visited, MapNode currNode, GeographicPoint goal,
                           HashMap<MapNode, MapNode> parentMap, PriorityQueue<MapNode> pq){
        //now we make our edge list
        Set<MapEdge> currNodeEdges = currNode.getEdges();
        for (MapEdge edge : currNodeEdges) {
            //get the neighbor at other end of edge
            MapNode aNeighbor = edge.getOtherNode(currNode);
            //if the neighbor has not bee visited, calculate distances
            if (!visited.contains(aNeighbor)){
                //this is where A* may be different
                //update aNeighbor's distance:
                //with the currNode's weighted distance, and the edge's length to the neighbor being looked at
                double neighborDistFromCurr = currNode.getActualDistance() + edge.getLength() ;
                double neighborToGoal = aNeighbor.getLocation().distance(goal);
                //distance from aNeighbor to the goal, this is give by the .distance() method
                //the formal distance + the heuristic distance calculated by neighborToGoal
                //if this distance is less than the neighbor's, we update aNeighbor's distance, and add it to parent map
                //so we can reach it, and add aNeighbor to the queue
                if (neighborToGoal+neighborDistFromCurr < aNeighbor.getDistance()) {
                    aNeighbor.setDistance(neighborToGoal+neighborDistFromCurr);
                    aNeighbor.setActualDistance(neighborDistFromCurr);
                    parentMap.put(aNeighbor, currNode);
                    pq.add(aNeighbor); //arbitrary added to queue data structure
                }
            }
        }
        //for each edge connected to the currNode
    }

    public List<GeographicPoint> tspPath(HashSet<GeographicPoint> tspNodeSet,
                                         GeographicPoint start) {
        //to track comparable shortest distances between edges
        LinkedList<GeographicPoint> optimalPath = new LinkedList<>();
        HashSet<MapNode> visited = new HashSet<>();
        MapNode closestNextNode = null;
        Double currMinDistance = Double.POSITIVE_INFINITY;

        MapNode currNode = pointNodeMap.get(start); //make currNode begin with start
        for (int i = 0; i < tspNodeSet.size()+1; i++) {
            for (GeographicPoint gp : tspNodeSet){
                MapNode nextNode = pointNodeMap.get(gp);
                //get an endNode from currNode according to current edge we are looking at
                //get the distance of the edge using adjusted aStarSearch Aglorithm (serves no purpose for simplemapgraph)
                double edgeLength = aStarSearchExtension(currNode.getLocation(), nextNode.getLocation());
                //if the edge is less than current min and the edge doesn't lead to a node we already seen, do the following
                if (edgeLength < currMinDistance && !visited.contains(nextNode)) {
                    //set minEdge = aEdge (this always happens atleast 1 time because currMinD is set to infinity)
                    closestNextNode = nextNode;
                    //set currMinD to the edgelength
                    currMinDistance = edgeLength;
                  //  tspNodeSet.remove(gp); //so that we don't have to
                }
            }
            currMinDistance = Double.POSITIVE_INFINITY; //reset our currMinDistance between nodes
            optimalPath.add(currNode.getLocation());
            currNode = closestNextNode; //the currNode is now the node at the end of our edge
            visited.add(currNode); //we can add the currNode to visited now
        }
        //so I'm getting smarter, but where is the heuristic?
        optimalPath.addLast(start);
        //because the last node we visit is going to be directly to the start again.
        return optimalPath;
    }

    //refactored out
    private double reconstructPathExtension(HashMap<MapNode,
            MapNode> parentMap, MapNode start, MapNode goal) {
        //we only need distance now
        //LinkedList<GeographicPoint> path = new LinkedList<GeographicPoint>();
        MapNode current = goal;
        MapNode next = null;
        double distance= 0;
        while (!current.equals(start)) {
            //we are not needing this path, we only need distance now
            //path.addFirst(current.getLocation());
            next = parentMap.get(current);
            //ystem.out.println(current.getLocation().distance(next.getLocation()));
            //track distance of the current nodes
            distance = distance + current.getLocation().distance(next.getLocation());
            current = next;
        }
        // we only need distance now
        //path.addFirst(start.getLocation());
        return distance;
    }
    //refactored out
    public boolean checkPoints(GeographicPoint start, GeographicPoint goal) {
        // Setup - check validity of inputs
        if (start == null || goal == null)
            throw new NullPointerException("Cannot find route from or to null node");
        MapNode startNode = pointNodeMap.get(start);
        MapNode endNode = pointNodeMap.get(goal);
        if (startNode == null) {
            System.err.println("Start node " + start + " does not exist");
            return false;
        }

        if (endNode == null) {
            System.err.println("End node " + goal + " does not exist");
            return false;
        }
        return true;
    }


    public static void main(String[] args)
    {
        MapGraph simpleTestMap = new MapGraph();
        GraphLoader.loadRoadMap("data/testdata/simpletest.map", simpleTestMap);

        GeographicPoint testStart = new GeographicPoint(4.0, 1.0); //this is startnode in tsp function
        GeographicPoint testPoint1 = new GeographicPoint(5.0, 1.0);
        GeographicPoint testPoint2= new GeographicPoint(4.0, 0.0);
        GeographicPoint testPoint3 = new GeographicPoint(5.0, 0.0);

        HashSet<GeographicPoint>tsp  = new HashSet<GeographicPoint>(){{
           add(testPoint1);
           add(testPoint2);
           add(testPoint3);
        }};

        System.out.println(simpleTestMap.tspPath(tsp,testStart));   //so far, print statements inside function


    }

}