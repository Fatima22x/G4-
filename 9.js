// CPCS 324 Algorithms & Data Structures 2
// Graph data structure starter - Transitive Closure Package
// 2018, Dr. Muhammad Al-Hashimi

// -----------------------------------------------------------------------
// simple graph object with linked-list edge implementation and minimal fields
// extra vertex and edge member fields and methods to be added later as needed
//


var _v = [], _e = [],       // exercise8_4_1.js
    _v2 = [], _e2 = [],     // exercise8_4_7.js
    _v3 = [], _v3 = [];     // ksa12cities.js

    
// -----------------------------------------------------------------------
function _main()
{
    // First graph for exercise8_4_1
    var g = new Graph();
    g.label = "Exercise 8.4: 1 (Levitin, 3rd edition)";
    g.digraph = true;
    g.readGraph(_v, _e);
    g.printGraph();

    g.connectedComp = g.topoSearch(g.dfs);
    document.write("<p>dfs_push: ", g.dfs_push, "</p>");

    // report connectivity status if available
    document.write("<p> ", g.componentInfo(), "</p>");

    g.makeAdjMatrix();

    // Output DFS-based TC matrix
    g.dfsTC();
    document.write("<p>TC matrix by DFS:</p>");
    for (var i = 0; i < g.DfsTC.length; i++)
    {
        document.write(g.DfsTC[i], "</br>");

    }

    // Output Warshall matrix
    g.warshallFloyd();
    document.write("<p>TC matrix by Warshall-Floyd:</p>");
    if (g.adjMatrix[0][0] != undefined)
    {
        for (var i = 0; i < g.adjMatrix.length; i++)
        {
            for (var j = 0; j < g.adjMatrix.length; j++)
            {
                document.write(g.adjMatrix[i][j].tc, j == g.adjMatrix.length - 1 ? " " : ",");
            }
            document.write("</br>");
        }
    }
    else 
    {
        document.write("Graph is not DAG, TC does not apply.")
    }  
    // Report if DAG or not
    document.write("<p>DAG: ", g.isDAG(), "</br>");



    // Second graph for exercise8_4_7
    var g2 = new Graph();
    g2.label = "Exercise 8.4: 1 (Levitin, 3rd edition)";
    g2.digraph = true;
    g2.readGraph(_v2, _e2);

    // Output Warshall matrix
    g2.warshallFloyd();
    document.write("<p>TC matrix by Warshall-Floyd Exercise 8.4: 7</p>");
    if (g2.adjMatrix[0][0].tc != undefined)
    {
        for (var i = 0; i < g2.adjMatrix.length; i++)
        {
            for (var j = 0; j < g2.adjMatrix.length; j++)
            {
                document.write(g2.adjMatrix[i][j].tc, j == g2.adjMatrix.length - 1 ? " " : ",");
            }
            document.write("</br>");
        }
    } 
    else 
    {
        document.write("Graph is not DAG, TC does not apply.")
    }  

    // Output FLoyd matrix
    if (g2.adjMatrix[0][0].dist != undefined)
    {
        document.write("<p>Distance matrix Exercise 8.4: 7</p>");
        for (var i = 0; i < g2.adjMatrix.length; i++)
        {
            for (var j = 0; j < g2.adjMatrix.length; j++)
            {
                document.write(g2.adjMatrix[i][j].dist, j == g2.adjMatrix.length - 1 ? " " : ",");
            }
            document.write("</br>");
        }
    }
    else 
    {
        document.write("Graph is not Weighted.")
    }  

    // Report if DAG or not
    document.write("<p>DAG Exercise 8.4: 7 </br>", g2.isDAG(), "</br>");



    // Third graph for ksa12cities
    var g3 = new Graph();
    g3.label = "KSA 12 Cities Dataset";
    g3.digraph = false;
    g3.readGraph(_v3, _e3);

    if (g3.weighted && !g3.digraph)
    {
        g3.prim();

        // output the Prim matrix
        document.write("<p>MST for " + g3.label + "</p>");
        var path = " " + g3.Prim[0].v.label + " ";
        var total = 0;
        for (var i = 0; i < g3.Prim.length; i++)
        {
            path = path + " " + g3.Prim[i].u.label + ", ";
            document.write("The path weight from " + g3.Prim[i].v.label + " to " + g3.Prim[i].u.label + " is: " +
                g3.Prim[i].w, " </br>");
            total = total + g3.Prim[i].w;
        }

        document.write("</br>Shortest path: " + path + " </br>Total weight: " + total);

    }

}




// -----------------------------------------------------------------------

function Vertex(v)
{
    // published docs section (ref. assignment page)
    // for this section, strip line comments
    // no JSDOC comments in this section

    // base base property fields

    this.label = v.label;
    this.visit = false;
    this.adjacent = new List();

    // base member methods

    this.adjacentById = adjacentByIdImpl;

    // --------------------
    // more student fields next


    // --------------------
    // more student methods next

    this.insertAdjacent = insertAdjacentImpl;
    this.vertexInfo = vertexInfoImpl;
    this.incidentEdges = incidentEdgesImpl;
}

// -----------------------------------------------------------------------

function Edge(vert_i, weight)
{
    // published docs section (ref. assignment page)
    // for this section, strip line comments
    // no JSDOC comments in this section


    // base property fields

    this.target_v = vert_i;
    if (weight === undefined)
    {
        this.weight = null;
    }
    else
    {
        this.weight = weight;
    }

    // base member methods


    // --------------------
    // more student fields next


    // --------------------
    // more student methods next

}


// -----------------------------------------------------------------------

function Graph()
{
    // published docs section (ref. assignment page)
    // for this section, strip line comments
    // no JSDOC comments in this section


    // base property fields

    this.vert = [];
    this.nv = 0;
    this.ne = 0;
    this.digraph = true;
    this.dfs_push = [];
    this.bfs_order = [];
    this.label = "";
    this.weighted = false;
    this.connectedComp = 0;
    this.adjMatrix = [];


    // base member methods

    this.readGraph = better_input;
    this.printGraph = printGraphImpl;
    this.listVerts = listVertsImpl;


    // --------------------
    // more student fields next

    this.DfsTC = [];
    // this.Warshall = [];
    // this.Floyd = [];
    this.Prim = [];

    // --------------------
    // more student methods next 

    this.addEdge = addEdgeImpl3;
    this.dfs = dfsImpl;
    this.bfs = bfsImpl;
    this.topoSearch = topoSearchImpl;
    this.makeAdjMatrix = makeAdjMatrixImpl3;
    this.isConnected = isConnectedImpl;
    this.componentInfo = componentInfoImpl;
    this.makeGraph = makeGraphImpl;

    // transitive closure package (requirements in line comments, to be removed and replaced by JSDOCs) 

    this.hasPath = hasPathImpl;
    this.shortestPath = shortestPathImpl;
    this.isDAG = isDAGImpl;
    this.warshallFloyd = warshallFloydImpl2;
    this.dfsTC = dfsTCImpl;
    this.prim = primImpl;

}


// -----------------------------------------------------------------------
// functions used by methods of Graph and ancillary objects

// -----------------------------------------------------------------------
// begin student code section
// -----------------------------------------------------------------------

// transitive closure package 

/**
 * @description // Checks if there's a path between vertices v_i, v_j
 * @param {integer} v_i source vertex id
 * @param {integer} v_j target vertex id
 * @author Fatima Falath
 * @returns {boolean} true if path exists between vertices v_i, v_j
 */

function hasPathImpl(v_i, v_j)
{
    return (this.DfsTC[v_i][v_j] === 1 ? true : false);
}

/**
 * @description Finds the shortest directed or udirected edge between vertices
 * @method dfsTC
 * @author Fatima Falath
 * @returns DfsTC matrix 
 */

function dfsTCImpl()
{
    for (var i = 0; i < this.nv; i++)
    {
        var v = this.vert[i];

        for (var j = 0; j < this.nv; j++)
        {
            this.vert[j].visit = false;
        }

        this.DfsTC[i] = [];
        for (var j = 0; j < this.nv; j++)
        {
            this.DfsTC[i][j] = 0;
        }

        var w = v.adjacentById();
        for (var k = 0; k < w.length; k++)
        {
            this.dfs(w[k]);
        }

        for (var m = 0; m < this.nv; m++)
        {
            if (this.vert[m].visit)
            {
                this.DfsTC[i][m] = 1;
            }
        }
    }
}

/**
 * @description obselote, returns TC matrix and distance matrix for diagraph based on dfs
 * @method warshallFloyd
 * @author Fatima Falath
 * @returns {object[], object[]} Warshall matrix, Floyd matrix
 */

function warshallFloydImpl()
{
    this.makeAdjMatrix();
    for (var i = 0; i < this.adjMatrix.length; i++)
    {
        this.Floyd[i] = this.adjMatrix[i].slice();
        this.Warshall[i] = this.adjMatrix[i].slice();

        // check if there is relation between vertices  
        for (var j = 0; j < this.nv; j++)
        {
            if (this.adjMatrix[i][j] == 0 && (i != j))
            {
                this.Floyd[i][j] = Infinity;
            }
        }
    }

    for (var k = 0; k < this.Warshall.length; k++)
    {
        for (var i = 0; i < this.Warshall.length; i++)
        {
            for (var j = 0; j < this.Warshall.length; j++)
            {
                this.Warshall[i][j] = (this.Warshall[i][j] ||
                    (this.Warshall[i][k] && this.Warshall[k][j])) ? 1 : 0;

                this.Floyd[i][j] =
                    Math.min(this.shortestPath(i, j), (this.shortestPath(i, k) + this.shortestPath(k, j)));
            }
        }

    }

    for (var i = 0; i < this.Floyd.length; i++)
    {
        for (var j = 0; j < this.Floyd.length; j++)
        {
            if (this.Floyd[i][j] == Infinity)
            {
                this.Floyd[i][j] = 0;
            }
        }
    }


}

/**
 * @description inserts .tc field in adjacency matrix if digraph, and .dist if weighted
 * @method warshallFloyd2
 * @author Fatima Falath
 * @returns {object[]} Adjacency matrix with new fields
 */

function warshallFloydImpl2()
{
    this.makeAdjMatrix();
    var newAdjMatrix = [];

    for (var i = 0; i < this.adjMatrix.length; i++)
    {
        newAdjMatrix[i] = this.adjMatrix[i].slice();

        // check if there is relation between vertices  
        if (this.weighted)
        {
            for (var j = 0; j < this.nv; j++)
            {
                if (this.adjMatrix[i][j] == 0 && (i != j))
                {
                    newAdjMatrix[i][j] = {"dist": Infinity};
                }
            }
        }
        
    }

    for (var k = 0; k < this.adjMatrix.length; k++)
    {
        for (var i = 0; i < this.adjMatrix.length; i++)
        {
            for (var j = 0; j < this.adjMatrix.length; j++)
            {
                this.adjMatrix[i][j] = {"tc": this.adjMatrix[i][j]};
                this.adjMatrix[i][j].tc = (this.adjMatrix[i][j].tc ||
                    (this.adjMatrix[i][k].tc && this.adjMatrix[k][j].tc)) ? 1 : 0;

                this.adjMatrix[i][j].dist =
                    Math.min(this.adjMatrix[i][j].dist, (this.adjMatrix[i][k].dist + this.adjMatrix[k][j].dist));
                
                if (this.weighted) 
                {
                    newAdjMatrix[i][j] = {
                            "dist": this.adjMatrix[i][j].dist
                        }
                }

                if (this.isDAG())
                {
                    newAdjMatrix[i][j] = {
                        "tc": this.adjMatrix[i][j].tc,
                    }
                }
                
            }
        }

    }
    
    this.adjMatrix = newAdjMatrix.map(function(arr)
    {
        return arr.slice();
    });

}

/**
 * @description return distance of shortest path between v_i, v_j in weighted graph 
 * @method shortestPath
 * @param {integer} v_i  Source vertex id
 * @param {integer} v_j  Target vertex id
 * @author Fatima Falath
 * @returns {edge} shortest path between two vertices
 */

function shortestPathImpl(v_i, v_j)
{   
    return (this.adjMatrix[v_i][v_j].dist);
}

/**
 * @description determines if graph is Directed Acyclic Graph
 * @method isDAG 
 * @author Fatima Falath
 * @returns {boolean} true if acyclic digraph
 */

function isDAGImpl()
{
    this.dfsTC();
    for (var i = 0; i < this.DfsTC.length; i++)
    {
        if (this.hasPath(i, i))
        {
            return false;
        }
    }
    return true;
}

/**
 * @description Finds Minimum Spanning Tree of graph
 * @method prim 
 * @author Fatima Falath
 * @returns {object[]} array of json objects
 */

function primImpl()
{

    var tree = [];

    // mark vertices unvisited
    for (var l = 0; l < this.nv; l++)
    {
        this.vert[l].visit = false;
    }

    // initiate first value as visited
    tree[0] = this.vert[0];
    this.vert[0].visit = true;

    var min = Infinity; // to find any number less than infinity
    for (var i = 0; i < this.nv; i++)
    {

        for (var j = 0; j < tree.length; j++)
        {

            var edges = tree[j].incidentEdges();

            for (var k = 0; k < edges.length; k++)
            {
                if (!this.vert[edges[k].adjVert_i].visit && edges[k].edgeWeight < min)
                {
                    this.Prim[i] =
                        (
                        {
                            v: tree[j],
                            u: this.vert[edges[k].adjVert_i],
                            w: edges[k].edgeWeight
                        });
                    min = this.Prim[i].w;
                }
            }
        }

        var n = this.Prim.length;
        tree[tree.length] = this.Prim[n - 1].u;
        this.Prim[n - 1].u.visit = true;
        min = Infinity;


    }
}

// -----------------------------------------------------------------------
// published docs section (ref. assignment page)
// use starter6-based P1M1 code as-is (fixes/improvements OK)
// no JSDOC comments in this section (docs already published)
// -----------------------------------------------------------------------

// --------------------
function addEdgeImpl(u_i, v_i) // obsolete, replaced by addEdgeImpl3() below
{
    // fetch edge vertices using their id, where u: source vertex, v: target vertex
    var u, v;
    u = this.vert[u_i];
    v = this.vert[v_i];

    // insert (u,v), i.e., insert v (by id) in adjacency list of u
    u.adjacent.insert(v_i);

    // insert (v,u) if undirected graph (repeat above but reverse vertex order)
    if (!this.digraph)
    {
        v.adjacent.insert(u_i);
    }
}

// --------------------
function addEdgeImpl2(u_i, v_i, weight) // obsolete, replaced by addEdgeImpl3() below
{
    // fetch vertices using their id, where u: edge source vertex, v: target vertex
    var u, v;
    u = this.vert[u_i];
    v = this.vert[v_i];

    // insert (u,v), i.e., insert v in adjacency list of u
    // (first create edge object using v_i as target, then pass object)
    var e1, e2;
    e1 = new Edge(v_i, weight);
    u.adjacent.insert(e1);

    // insert (v,u) if undirected graph (repeat above but reverse vertex order)
    if (!this.digraph)
    {
        e2 = new Edge(u_i, weight);
        v.adjacent.insert(e2);
    }

}

// --------------------
function addEdgeImpl3(u_i, v_i, weight)
{
    // fetch vertices using their id, where u: edge source vertex, v: target vertex
    var u, v;
    u = this.vert[u_i];
    v = this.vert[v_i];

    // insert (u,v), i.e., insert v in adjacency list of u
    u.insertAdjacent(v_i, weight);

    // insert (v,u) if undirected graph (repeat above but reverse vertex order)
    if (!this.digraph)
    {
        v.insertAdjacent(u_i, weight);
    }

}

// --------------------
function adjacentByIdImpl()
{
    var out = [];
    var e = this.adjacent.traverse();
    for (var i = 0; i < e.length; i++)
    {
        out[i] = e[i].target_v;
    }
    return out;
}

// --------------------
function better_input(v, e)
{
    // set vertex and edge count fields
    this.nv = v.length;
    this.ne = e.length;


    // input vertices into internal vertex array
    for (i = 0; i < this.nv; i++)
    {
        this.vert[i] = new Vertex(v[i]);
    }

    // input vertex pairs from edge list input array
    // remember to pass vertex ids to add_edge()
    for (i = 0; i < this.ne; i++)
    {
        if (e[i].w != null)
        {
            this.addEdge(e[i].u, e[i].v, e[i].w);
            continue;
        }
        this.addEdge(e[i].u, e[i].v);
    }

    // double edge count if graph undirected
    if (!this.digraph)
    {
        this.ne = 2 * e.length;
    }

    // check if the graph is weighted or not 
    if (!(e[0].w === undefined))
    {
        this.weighted = true;
    }


}

// --------------------
function better_output()
{
    document.write("<p>GRAPH {", this.label, "} ", this.weighted ? "" : "UN", "WEIGHTED, ",
        this.digraph ? "" : "UN", "DIRECTED - ", this.nv, " VERTICES, ",
        this.ne, " EDGES:</p>");

    // list vertices
    this.list_vert();

    // report connectivity status if available
    connected = (this.connectedComp === 0 ? true : false);
    document.write("<p>", connected ? "no connectivity info" : ("DISCONNECTED " + this.connectedComp), "</p>");


}

// --------------------
function bfsImpl(v_i)
{
    var i;

    // get vertex v by its id
    var v = this.vert[v_i];

    // process v 
    v.visit = true;
    this.bfs_order[this.bfs_order.length] = v_i;

    // initialize queue with v
    var queue = new Queue();
    queue.enqueue(v);

    // while queue not empty
    while (!queue.isEmpty())
    {
        // dequeue and process a vertex, u
        var u = queue.dequeue();

        // queue all unvisited vertices adjacent to u
        var w = u.adjacentById();
        for (i = 0; i < w.length; i++)
        {
            if (!this.vert[w[i]].visit)
            {
                {
                    this.vert[w[i]].visit = true;
                    queue.enqueue(this.vert[w[i]]);
                    this.bfs_order[this.bfs_order.length] = w[i];
                }
            }
        }
    }
}

// --------------------
function componentInfoImpl()
{
    // report connectivity status if available
    if (this.isConnected() === true)
    {
        return "no connectivity info";
    }
    return this.connectedComp === 1 ? "CONNECTED" : ("DISCONNECTED " + this.connectedComp);
}

// --------------------
function dfsImpl(v_i)
{
    // get landing vert by id then process
    var v = this.vert[v_i];
    v.visit = true;
    this.dfs_push[this.dfs_push.length] = v_i; //insert into dfs_push array

    // recursively traverse unvisited adjacent vertices
    var w = v.adjacentById();
    for (var i = 0; i < w.length; i++)
    {
        if (!this.vert[w[i]].visit)
        {
            this.dfs(w[i]);
        }

    }
}

// --------------------
function incidentEdgesImpl()
{
    var edgeNodes = [];
    var e = this.adjacent.traverse();
    for (var i = 0; i < e.length; i++)
    {
        edgeNodes[i] = {
            "adjVert_i": e[i].target_v,
            "edgeLabel": "",
            "edgeWeight": e[i].weight
        };
    }
    return edgeNodes;
}

// --------------------
function insertAdjacentImpl(v_i, weight)
{
    var e = new Edge(v_i);

    if (!(weight === undefined))
    {
        e.weight = weight;
    }

    this.adjacent.insert(e);
}

// --------------------
function isConnectedImpl()
{
    return (this.connectedComp === 0 ? true : false);
}

// --------------------
function listVertsImpl()
{
    var i, v; // local vars
    for (i = 0; i < this.nv; i++)
    {
        v = this.vert[i];
        document.write(v.vertexInfo(i));
    }
}

// --------------------
function makeAdjMatrixImpl() //obselote, replaced by makeAdjMatrixImpl3 below
{
    // initially zero the adjecency matrix
    for (var i = 0; i < this.nv; i++)
    {
        this.adjMatrix[i] = [];
        for (var j = 0; j < this.nv; j++)
        {
            this.adjMatrix[i][j] = 0;
        }
    }

    // for each vertex, set 1 for each adjacency
    var v, w, m;
    for (var i = 0; i < this.nv; i++)
    {
        v = this.vert[i];
        w = v.adjacentById();
        for (var j = 0; j < w.length; j++)
        {
            this.adjMatrix[i][w[j]] = 1;
        }
    }

}

// --------------------
function makeAdjMatrixImpl2() //obselote, replaced by makeAdjMatrixImpl3 below
{
    // initially create row elements and zero the adjacency matrix
    for (var i = 0; i < this.nv; i++)
    {

        // get vertex
        var v = this.vert[i];

        // create and init the corresponding row
        this.adjMatrix[i] = [];
        for (var j = 0; j < this.nv; j++)
        {
            this.adjMatrix[i][j] = 0;
        }

        // process adjacent vertices: get edge node, set value for each
        var w = v.adjacent.traverse();
        for (j = 0; j < w.length; j++)
        {
            this.adjMatrix[i][w[j].target_v] = this.weighted ? w[j].weight : 1;
        }
    }
}

// --------------------
function makeAdjMatrixImpl3()
{
    // initially create row elements and zero the adjacency matrix
    for (var i = 0; i < this.nv; i++)
    {

        // get vertex
        var v = this.vert[i];

        // create and init the corresponding row
        this.adjMatrix[i] = [];
        for (var j = 0; j < this.nv; j++)
        {
            this.adjMatrix[i][j] = 0;
        }

        // process adjacent vertices: get edge node, set value for each
        var w = v.incidentEdges();
        for (j = 0; j < w.length; j++)
        {
            this.adjMatrix[i][w[j].adjVert_i] = this.weighted ? w[j].edgeWeight : 1;
        }
    }
}


// To be implemented later
function makeGraphImpl(n, m, w)
{

}

// --------------------
function topoSearchImpl(fun)
{
    var i, c = 0;

    // mark all vertices unvisited
    for (i = 0; i < this.nv; i++)
    {
        this.vert[i].visit = false;
    }

    // traverse unvisited connected component 	
    for (i = 0; i < this.nv; i++)
    {
        if (!this.vert[i].visit)
        {
            switch (fun)
            {
                case this.dfs:
                    (c++, this.dfs(i));
                    break;
                case this.bfs:
                    this.bfs(i);
                    break;
            }
        }
    }
    return c;
}

// --------------------
function printGraphImpl()
{
    document.write("<p>GRAPH {", this.label, "} ", this.weighted ? "" : "UN", "WEIGHTED, ",
        this.digraph ? "" : "UN", "DIRECTED - ", this.nv, " VERTICES, ",
        this.ne, " EDGES:</p>");

    // report connectivity status if available
    document.write("<p> ", this.componentInfo(), "</p>");

    // list vertices
    this.listVerts();




}

// --------------------
function vertexInfoImpl(vert_i)
{
    return ("VERTEX:" + vert_i + " {" + this.label + "} - VISIT: " + this.visit + " - ADJACENCY: " + this.adjacentById() + "<br>");
}
