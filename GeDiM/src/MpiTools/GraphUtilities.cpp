#include "GraphUtilities.hpp"

namespace Gedim
{
  // ***************************************************************************
  void GraphUtilities::DepthFirstSearch(const unsigned int& vertex,
                                        const std::vector<std::vector<unsigned int> >& graphAdjacency,
                                        std::vector<bool>& visited,
                                        std::list<unsigned int>& visitedVertices) const
  {
    // Mark the current node as visited and print it
    visited[vertex] = true;
    visitedVertices.push_back(vertex);

    // Recur for all the vertices adjacent to this vertex
    for (const unsigned int& adj_v : graphAdjacency[vertex])
    {
      if (visited[adj_v])
        continue;

      DepthFirstSearch(adj_v,
                       graphAdjacency,
                       visited,
                       visitedVertices);
    }
  }
  // ***************************************************************************
  std::vector<std::vector<unsigned int>> GraphUtilities::ComputeAdjacencyTranspose(const unsigned int& graphNumVertices,
                                                                                   const std::vector<std::vector<unsigned int> >& graphAdjacency) const
  {

    std::vector<std::vector<unsigned int>> graphT_Adjacency(graphNumVertices);

    std::vector<std::list<unsigned int>> adj(graphNumVertices);
    for (unsigned int v = 0; v < graphNumVertices; v++)
    {
      for (const unsigned int& adj_v : graphAdjacency[v])
        adj[adj_v].push_back(v);
    }

    for (unsigned int v = 0; v < graphNumVertices; v++)
      graphT_Adjacency[v] = std::vector<unsigned int>(adj[v].begin(),
                                                      adj[v].end());

    return graphT_Adjacency;
  }
  // ***************************************************************************
  void GraphUtilities::FillOrder(const unsigned int& v,
                                 const std::vector<std::vector<unsigned int> >& graphAdjacency,
                                 std::vector<bool>& visited,
                                 std::stack<int>& stack) const
  {
    // Mark the current node as visited and print it
    visited[v] = true;

    // Recur for all the vertices adjacent to this vertex
    for (const unsigned int& adj_v : graphAdjacency[v])
    {
      if (visited[adj_v])
        continue;

      FillOrder(adj_v,
                graphAdjacency,
                visited,
                stack);
    }

    // All vertices reachable from v are processed by now, push v
    stack.push(v);
  }
  // ***************************************************************************
  std::vector<std::vector<unsigned int>> GraphUtilities::GraphConnectivityToGraphAdjacency(const unsigned int& graphNumVertices,
                                                                                           const Eigen::MatrixXi& graphConnectivity) const
  {
    std::vector<std::vector<unsigned int>> adjacency(graphNumVertices);

    std::vector<std::list<unsigned int>> adj(graphNumVertices);
    for (unsigned int e = 0; e < graphConnectivity.cols(); e++)
      adj[graphConnectivity(0, e)].push_back(graphConnectivity(1, e));

    for (unsigned int v = 0; v < graphNumVertices; v++)
      adjacency[v] = std::vector<unsigned int>(adj[v].begin(),
                                               adj[v].end());

    return adjacency;
  }
  // ***************************************************************************
  Eigen::MatrixXi GraphUtilities::GraphAdjacencyToGraphConnectivity(const unsigned int& graphNumEdges,
                                                                    const std::vector<std::vector<unsigned int>>& graphAdjacency) const
  {
    Eigen::MatrixXi graphConnectivity(2, graphNumEdges);

    unsigned int e = 0;
    for (unsigned int v = 0; v < graphAdjacency.size(); v++)
    {
      for (const unsigned int& adj_v : graphAdjacency[v])
        graphConnectivity.col(e++)<< v, adj_v;
    }

    return graphConnectivity;
  }
  // ***************************************************************************
  std::vector<std::vector<unsigned int>> GraphUtilities::ComputeStronglyConnectedComponents(const unsigned int& graphNumVertices,
                                                                                            const std::vector<std::vector<unsigned int> >& graphAdjacency) const
  {
    std::stack<int> stack;

    // Mark all the vertices as not visited (For first DFS)
    std::vector<bool> visited(graphNumVertices, false);

    // Fill vertices in stack according to their finishing times
    for (unsigned int i = 0; i < graphNumVertices; i++)
    {
      if (visited[i])
        continue;

      FillOrder(i,
                graphAdjacency,
                visited,
                stack);
    }

    // Create a reversed graph
    const std::vector<std::vector<unsigned int>> graphT_Adjacency = ComputeAdjacencyTranspose(graphNumVertices,
                                                                                              graphAdjacency);

    // Mark all the vertices as not visited (For second DFS)
    for(unsigned int i = 0; i < graphNumVertices; i++)
      visited[i] = false;

    std::list<std::vector<unsigned int>> connectedComponents;

    // Now process all vertices in order defined by Stack
    while (stack.empty() == false)
    {
      // Pop a vertex from stack
      int v = stack.top();
      stack.pop();

      // Print Strongly connected component of the popped vertex
      if (visited[v])
        continue;

      std::list<unsigned int> connectedComponent;
      DepthFirstSearch(v,
                       graphT_Adjacency,
                       visited,
                       connectedComponent);

      connectedComponents.push_back(std::vector<unsigned int>(connectedComponent.begin(),
                                                              connectedComponent.end()));
    }

    return std::vector<std::vector<unsigned int>>(connectedComponents.begin(),
                                                  connectedComponents.end());
  }
  // ***************************************************************************
}