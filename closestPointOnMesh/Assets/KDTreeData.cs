//#define DEBUG_KD
using UnityEngine;
using System.Collections;
using System.Collections.Generic;
using UniRx;

public class KDTreeData {

    static Dictionary<int, Node> dict = new Dictionary<int, Node>();

    //[SerializeField]
    //bool drawMeshTreeOnStart=false;

    public class Node
    {
        public float partitionCoordinate;
        public int partitionAxis;

        public Node positiveChild;
        public Node negativeChild;

        public int[] triangles;

        public Bounds bounds;

        public Vector3 tempClosestPoint;
    };

    private int triangleCount;
    private int vertexCount;
    private Vector3[] vertices;
    private int[] tris;
    //private float[] areas;
    private Vector3[] triangleCentroid;
    private float[] triangleRadius;
    // private Vector3[] triangleNormals;
    // private int[] randomTris;
    // private float[] tempWeights;

    private bool[] needToCheck;
    private int[] triangleToClear;
    // private float totalArea = 0f;

    private Mesh mesh;

    public Node rootNode;

    public int ID { get; private set; }

    public KDTreeData(Mesh mesh) {

        Processing = true;

        this.mesh = mesh;

        ID       = mesh.GetInstanceID();

        tris     = mesh.triangles;
        vertices = mesh.vertices;
    }

    public bool Processing = false;

    public void Build() {

        vertexCount = vertices.Length;
        triangleCount = tris.Length / 3;

        //triangleNormals   = new Vector3 [triangleCount];
        //areas             = new float   [triangleCount];
        triangleCentroid    = new Vector3 [triangleCount];
        triangleRadius      = new float   [triangleCount];
        //tempWeights       = new float   [triangleCount];

        needToCheck         = new bool [triangleCount];
        triangleToClear     = new int  [triangleCount];

        for (int i = 0; i < tris.Length; i += 3) {

            int index_a = tris[i + 0];
            int index_b = tris[i + 1];
            int index_c = tris[i + 2];

            Vector3 a = vertices[index_a];
            Vector3 b = vertices[index_b];
            Vector3 c = vertices[index_c];

            Vector3 ab = b - a;
            Vector3 ac = c - a;

            //Vector3 normal = Vector3.Cross((ab).normalized, (ac).normalized).normalized;

            //float area = Mathf.Abs(Vector3.Cross(ab, ac).magnitude / 2f);
            //areas[i / 3] = area;
            // totalArea += area;

            //triangleNormals[i / 3] = normal;

            Vector3 mean = (a + b + c) / 3f;
            triangleCentroid[i / 3] = mean;

            triangleRadius[i / 3] =
                Mathf.Max(
                    Vector3.Distance(a, mean),
                    Vector3.Distance(b, mean),
                    Vector3.Distance(c, mean)
                ) + float.Epsilon + 0.00000001f;

            needToCheck[i / 3] = true;
        }

        // computationMatrices = new Matrix[triangleCount+1];

        /*
        for (int i=0; i < triangleCount; i++)
        {
            computationMatrices[i] = new Matrix(i, 3, 0d);
        }*/

        //if (!drawMeshTreeOnStart)
        BuildTriangleTree();

        listPool.Clear();

        Processing = true;
    }

    Stack<Node> nodesToProcess = new Stack<Node>(1000);
    
    // bruteforce
    /*public Vector3 ClosestPointOnBruteForce(Vector3 to)
    {
        float currentMin = Mathf.Infinity;
        int currTriangle = 0;

        for (int i=0; i < triangleCount; i++)
        {
            
            //approx test
            if (MinDistanceToTriangleApprox(i, to) <= currentMin)
            {
                float newDist = ClosestDistanceOnTriangleSingle(i, to); //first test cache triangle

                if (newDist < currentMin)
                {
                    currentMin = newDist;
                    currTriangle = i;
                }
            }
        }

        Vector3 result;

        ClosestPointOnTriangleToPoint(
            ref vertices[tris[currTriangle]],
            ref vertices[tris[currTriangle + 1]],
            ref vertices[tris[currTriangle + 2]],
            ref to,
            out result);

        return result;
    }*/
    
    // List<int> tempTriangles = new List<int>();
    /*public Vector3 ClosestPointOnOptimizedNoRange(Vector3 to)
    {
        Vector3 vec = Vector3.one;

        ClosestPointOnOptimized(to, out vec, -1f);

        return vec;
    }

    public Vector3 ClosestPointOnOptimizedCached(Vector3 to)
    {
        Vector3 vec = Vector3.one;

        ClosestPointOnOptimized(to, out vec, -1f);

        return vec;
    }*/

    public int closestCachedTriangle = -1;
    public bool ClosestPointOnOptimized(Vector3 to, out Vector3 closestPoint, /*ref int triangleCacheIndex,*/ Transform transform, float range = -1f) 
    {
        
        to = transform.InverseTransformPoint(to);

        //tempTriangles.Clear();
        
        // clear queue
        nodesToProcess.Clear();
        
        float currentMin = 0;
        int currClosestTriangle = 0;
        int triangleToClearCount = 0;
        //maximum range is already defined
        if (range != -1f)
        {
            Vector3 scl = transform.lossyScale;
            range /= Mathf.Min(scl.x, Mathf.Min(scl.y, scl.z));

            currentMin = range;
        }
        /*else if(closestCachedTriangle != -1) //TODO TRIANGLES
        {

            int newClosestTriangle = node.triangles[i];

            int dividedIndex = newClosestTriangle / 3;

            if (needToCheck[dividedIndex])
            {
                needToCheck[dividedIndex] = false;
                triangleToClear[triangleToClearCount] = dividedIndex;
                triangleToClearCount++;

            currentMin = ClosestDistanceOnTriangleSingle(closestCachedTriangle, to); //first test cache triangle
        }*/
        else
        {
            currentMin = 9999;
        }

        /*if (range == -1f)
        {
            // Check random triangles for approximate search radius
            for (int i=0; i < randomTris.Length; i++)
            {
                int newClosestTriangle = randomTris[i];

                int dividedIndex = newClosestTriangle / 3;

                if (needToCheck[dividedIndex])
                {
                    needToCheck[dividedIndex] = false;
                    triangleToClear[triangleToClearCount] = dividedIndex;
                    triangleToClearCount++;
                }
                else
                    continue; //already checked, move along

                // approx test
                if (MinDistanceToTriangleApprox(newClosestTriangle, to) <= currentMin)
                {
                    // actual test
                    float newMin = ClosestDistanceOnTriangleSingle(newClosestTriangle, to);

                    // je ta trikotnik bližje?
                    if (newMin <= currentMin)
                    {
                        currentMin = newMin;
                        currClosestTriangle = newClosestTriangle;
                    }
                }
            }
        }*/

        rootNode.tempClosestPoint = rootNode.bounds.ClosestPoint(to);

        // push root tree
        nodesToProcess.Push(rootNode);

        int checkCount = 0;
        // tree.tempMinMax = PointDistanceFromPlane(tree.partitionPoint, tree.partitionNormal, to);
        // BSP search with pruning (don't visit areas which distance is more away than curr min distance
        while (nodesToProcess.Count > 0)
        {
            var node = nodesToProcess.Pop();

#if DEBUG_KD
            DebugDraw.DrawBounds(node.bounds, Color.red, 0f, transform);
#endif
            //if (Vector3.SqrMagnitude(node.tempClosestPoint - to) > currentMin * currentMin)
            //  continue;

            // pruning!
            if (node.triangles == null)
            {   
                // var   distToPlaneUnsigned = PointDistanceFromPlaneNoAbs(node.partitionCoordinate, node.partitionAxis, to);
                // float distAbs = Mathf.Abs(distToPlaneUnsigned);
                
                int   partitionAxis = node.partitionAxis;
                float partitionCoord = node.partitionCoordinate;

                Vector3 tempClosestPoint = node.tempClosestPoint;
                
                // Calculating closest distance to bounding box

                // Debug.Log("Distance:" + (tempClosestPoint[partitionAxis] - partitionCoord));
                
                //inside positive side, project on negative
                if ((tempClosestPoint[partitionAxis] - partitionCoord) >= 0)
                {
                    node.positiveChild.tempClosestPoint = tempClosestPoint;

                    tempClosestPoint[partitionAxis] = partitionCoord;
                    node.negativeChild.tempClosestPoint = tempClosestPoint;

                    //we are inside positive bound, we don't need to test for distance
                    nodesToProcess.Push(node.positiveChild);

                    if (Vector3.SqrMagnitude(node.negativeChild.tempClosestPoint - to) <= currentMin * currentMin
                    && (node.negativeChild.triangles != null && node.negativeChild.triangles.Length != 0))
                        nodesToProcess.Push(node.negativeChild);
                }
                else //inside negative side, project on positive
                {
                    node.negativeChild.tempClosestPoint = tempClosestPoint;

                    tempClosestPoint[partitionAxis] = partitionCoord;
                    node.positiveChild.tempClosestPoint = tempClosestPoint;

                    if (Vector3.SqrMagnitude(node.positiveChild.tempClosestPoint - to) <= currentMin * currentMin
                    && (node.positiveChild.triangles != null && node.positiveChild.triangles.Length != 0))
                    
                        nodesToProcess.Push(node.positiveChild);

                    //we are inside negative bound, we don't need to test for distance
                    nodesToProcess.Push(node.negativeChild);
                }

                ///Debug.Log("P:" + Vector3.SqrMagnitude(node.positiveChild.tempClosestPoint - to));
                ///Debug.Log("N:" + Vector3.SqrMagnitude(node.negativeChild.tempClosestPoint - to));

            }
            else
            {
                // preglejmo vse trikotnike, vzemimo najbližjega
                for (int i = 0; i < node.triangles.Length; i++)
                {
                    checkCount++;

                    int newClosestTriangle = node.triangles[i];

                    int dividedIndex = newClosestTriangle / 3;

                    if (needToCheck[dividedIndex])
                    {
                        needToCheck[dividedIndex] = false;
                        triangleToClear[triangleToClearCount] = dividedIndex;
                        triangleToClearCount++;
                    }
                    else
                        continue; //already checked, move along

                    // approx test
                    if (MinDistanceToTriangleApprox(newClosestTriangle, to) <= currentMin)
                    {

                        // actual test
                        float newMin = ClosestDistanceOnTriangleSingle(newClosestTriangle, to);
#if DEBUG_KD
                        DebugDraw.DrawTriangle(
                            transform.TransformPoint(vertices[tris[newClosestTriangle    ]]),
                            transform.TransformPoint(vertices[tris[newClosestTriangle + 1]]),
                            transform.TransformPoint(vertices[tris[newClosestTriangle + 2]]),
                            Color.yellow
                        );
#endif
                        // je ta trikotnik bližje?
                        if (newMin <= currentMin)
                        {
                            currentMin = newMin;
                            currClosestTriangle = newClosestTriangle;
                        }
                    }
                }
            }
        }

#if DEBUG_KD
        Debug.Log("checkCount:" + checkCount);
        Debug.Log("actualTrianglesTested:" + triangleToClearCount);
#endif
        while(triangleToClearCount > 0 ) 
        {
            triangleToClearCount--;
            needToCheck[triangleToClear[triangleToClearCount]] = true;            
        }

        if (checkCount == 0){ 
            closestPoint = Vector3.zero;
            return false;
        }

        /*foreach (bool b in needToCheck)
        {
            if (!b)
                Debug.Log("GOTYOU");
        }
        */
        Vector3 closest = Vector3.zero;

        ClosestPointOnTriangleToPoint(
            ref vertices[tris[currClosestTriangle]],
            ref vertices[tris[currClosestTriangle + 1]],
            ref vertices[tris[currClosestTriangle + 2]],
            ref to, 
            out closest);

#if DEBUG_KD
        DebugDraw.DrawTriangle(
            transform.TransformPoint(vertices[tris[currClosestTriangle    ]]), 
            transform.TransformPoint(vertices[tris[currClosestTriangle + 1]]), 
            transform.TransformPoint(vertices[tris[currClosestTriangle + 2]]),
            Color.red
            );
#endif

        closestPoint = transform.TransformPoint(closest);


        return true;
    }

    // temp
    Vector3 p1;
    Vector3 p2;
    Vector3 p3;
    Vector3 nearest;
    float ClosestDistanceOnTriangleSingle(int triangle, Vector3 to)
    {
        p1 = vertices[tris[triangle    ]];
        p2 = vertices[tris[triangle + 1]];
        p3 = vertices[tris[triangle + 2]];

        ClosestPointOnTriangleToPoint(ref p1, ref p2, ref p3, ref to, out nearest);

        return Vector3.Magnitude(to - nearest);
    }
    /*
    float ClosestDistanceOnTriangleSingleSqr(int triangle, Vector3 to)
    {
        p1 = vertices[tris[triangle]];
        p2 = vertices[tris[triangle + 1]];
        p3 = vertices[tris[triangle + 2]];

        ClosestPointOnTriangleToPoint(ref p1, ref p2, ref p3, ref to, out nearest);

        return (to - nearest).sqrMagnitude;
    }*/
    /*
    
    float Sqr(float x) {
        return x * x;
    }

    float ClosestDistanceToBox(Vector3 p, Bounds r) {
        float d = 0;
        if      (p.x<r.min.x) d += Sqr(r.min.x - p.x);
        else if (p.x>r.max.x) d += Sqr(p.x - r.max.x);
        
        if      (p.y<r.min.y) d += Sqr(r.min.y - p.y);
        else if (p.y>r.max.y) d += Sqr(p.y - r.max.y);
        
        if      (p.z<r.min.z) d += Sqr(r.min.z - p.z);
        else if (p.z>r.max.z) d += Sqr(p.z - r.max.z);
        
        return d;
    }*/

    /// <summary>
    /// Can be negative in weird cases, still you need to check
    /// </summary>
    /// <param name="triangle"></param>
    /// <param name="to"></param>
    /// <returns></returns>
    float MinDistanceToTriangleApprox(int triangle, Vector3 to)
    {
        int triangleIndex = triangle / 3;

        return Vector3.Distance(triangleCentroid[triangleIndex], to) - triangleRadius[triangleIndex];
    }

    Vector3 ClosestPointOnTriangle(int[] triangles, Vector3 to)
    {
        float shortestDistance = float.MaxValue;

        Vector3 shortestPoint = Vector3.zero;

        // Iterate through all triangles
        foreach (int triangle in triangles)
        {
            Vector3 p1 = vertices[tris[triangle + 0]];
            Vector3 p2 = vertices[tris[triangle + 1]];
            Vector3 p3 = vertices[tris[triangle + 2]];

            Vector3 nearest;

            ClosestPointOnTriangleToPoint(ref p1, ref p2, ref p3, ref to, out nearest);

            float distance = (to - nearest).sqrMagnitude;

            if (distance <= shortestDistance)
            {
                shortestDistance = distance;
                shortestPoint = nearest;
            }
        }

        return shortestPoint;
    }

    void BuildTriangleTree()
    {
        List<int> rootTriangles = new List<int>();

        for (int i = 0; i < tris.Length; i += 3)
        {
            rootTriangles.Add(i);
        }

        rootNode = new Node();
        rootNode.partitionAxis = -1;
        rootNode.bounds = MakeBounds(rootTriangles);

        
        RecursivePartition(rootTriangles, 0, rootNode, -1);
        

        // ChooseEvenlyRandomTriangles();
    }

    Bounds MakeBounds(List<int> triangles)
    {
        Vector3 maxExtents = new Vector3(-float.MaxValue, -float.MaxValue, -float.MaxValue);
        Vector3 minExtents = new Vector3( float.MaxValue,  float.MaxValue,  float.MaxValue);

        triangles.ForEach((triangle)=> {
        
            // partitionMeanPoint += centroids[triangle/3];

            minExtents.x = Mathf.Min(minExtents.x, Mathf.Min(vertices[tris[triangle]].x, Mathf.Min(vertices[tris[triangle + 1]].x, Mathf.Min(vertices[tris[triangle + 2]].x))));
            minExtents.y = Mathf.Min(minExtents.y, Mathf.Min(vertices[tris[triangle]].y, Mathf.Min(vertices[tris[triangle + 1]].y, Mathf.Min(vertices[tris[triangle + 2]].y))));
            minExtents.z = Mathf.Min(minExtents.z, Mathf.Min(vertices[tris[triangle]].z, Mathf.Min(vertices[tris[triangle + 1]].z, Mathf.Min(vertices[tris[triangle + 2]].z))));

            maxExtents.x = Mathf.Max(maxExtents.x, Mathf.Max(vertices[tris[triangle]].x, Mathf.Max(vertices[tris[triangle + 1]].x, Mathf.Max(vertices[tris[triangle + 2]].x))));
            maxExtents.y = Mathf.Max(maxExtents.y, Mathf.Max(vertices[tris[triangle]].y, Mathf.Max(vertices[tris[triangle + 1]].y, Mathf.Max(vertices[tris[triangle + 2]].y))));
            maxExtents.z = Mathf.Max(maxExtents.z, Mathf.Max(vertices[tris[triangle]].z, Mathf.Max(vertices[tris[triangle + 1]].z, Mathf.Max(vertices[tris[triangle + 2]].z))));
        });

        Bounds b = new Bounds();
        b.min = minExtents;
        b.max = maxExtents;
        return b;
    }

    Stack<List<int>> listPool = new Stack<List<int>>();

    List<int> GetList() {

        if (listPool.Count > 0)
        {
            var l = listPool.Pop();
            l.Clear();
            return l;
        }
        else
        {
            return new List<int>();
        }
    }

    void ReturnList(List<int> list)
    {
        listPool.Push(list);
    }

    void RecursivePartition(List<int> triangles, int depth, Node parent, int previousAxis)
    {
        int triCount = triangles.Count;

        float partitionCoordinate = 0f;
        
        // center of bounding box
        // Vector3 partitionMeanPoint = parent.bounds.center;
        Bounds parentBounds = parent.bounds;
        Vector3 parentBoundsSize = parentBounds.size;

        int partitionAxis = 0;
        float maxSize = parentBoundsSize.x;
        
        //is Y axis bigger?
        if (maxSize < parentBoundsSize.y) {
            partitionAxis = 1;
            maxSize = parentBoundsSize.y;
        }
        
        //is Z axis biggeR?
        if (maxSize < parentBoundsSize.z) {
            partitionAxis = 2;
            maxSize = parentBoundsSize.z;
        }

        /*
        // choosing axis, must alternate
        if (previousAxis != -1)
            partitionAxis = (previousAxis + 1) % 3;
        else
        if (extentsMagnitude.x >= extentsMagnitude.y && extentsMagnitude.x >= extentsMagnitude.z)
            partitionAxis = 0;
        else if (extentsMagnitude.y >= extentsMagnitude.x && extentsMagnitude.y >= extentsMagnitude.z)
            partitionAxis = 1;
        else
            partitionAxis = 2;
        */
        //partitionCoordinate = partitionMeanPoint[partitionAxis];

        //area that doesn't change when re-positioning splitting plane
        float sideArea = 2 * parentBoundsSize[(partitionAxis + 1) % 3] * parentBoundsSize[(partitionAxis + 2) % 3];

        //area that changes when re-positioning splitting plane
        float axisArea = 2 * parentBoundsSize[partitionAxis] * (parentBoundsSize[(partitionAxis + 1) % 3] + parentBoundsSize[(partitionAxis + 2) % 3]);

        float startPos = parentBounds.min[partitionAxis];
        float endPos   = parentBounds.max[partitionAxis];

        float minPosition = 0;
        float minHeuristic = System.Single.MaxValue;

        int areas = Mathf.Max(triangles.Count / 20, 5);
        
        float step = (endPos - startPos) / areas;

        // sweep test
        for (int i = 1; i < areas; i++) {

            float pos = startPos + (step * i);

            float h = Heuristic(triangles, pos, startPos, endPos, partitionAxis, axisArea, sideArea);

            if (h < minHeuristic) {

                minHeuristic = h;
                minPosition = pos;
            }
        }

        step /= 2f;

        // making partition, bisection
        // maximum 12 iterations!
        // works OK
        int k = 8;
        while (k > 0) {
            // desna razlika
            float right_h = Heuristic(triangles, minPosition + step, startPos, endPos, partitionAxis, axisArea, sideArea);

            // leva razlika 
            float left_h  = Heuristic(triangles, minPosition - step, startPos, endPos, partitionAxis, axisArea, sideArea);
            
            // Count(triangles, partitionCoordinate, partitionAxis, out posCount, out negCount);
            // int trianglesInBoth = triangles.Count - posCount + negCount;

            //moving right is positive, it makes heuristic smaller
            if (right_h < minHeuristic && right_h < left_h)
            {
                minPosition += step;
                minHeuristic = right_h;
            }
            //moving left is negative, it makes heuristic smaller
            else 
            if (left_h < minHeuristic && left_h < right_h)
            {
                minPosition -= step;
                minHeuristic = left_h;
            }
            else
                break;
            
            //minimize step by two
            step /= 2f;

            k--;
        }

        partitionCoordinate = minPosition;

        //partitionMeanPoint[partitionAxis] = partitionCoordinate;
        // DebugDraw.DrawMarker(transform.TransformPoint(partitionMeanPoint), 1f, Color.red, 10f, false);

        List<int> positiveTriangles = GetList();
        List<int> negativeTriangles = GetList();

        Split(triangles, partitionCoordinate, partitionAxis, positiveTriangles, negativeTriangles);

        ReturnList(triangles);

        parent.partitionAxis       = partitionAxis;
        parent.partitionCoordinate = partitionCoordinate;

        //Vector3 planeNormal = Vector3.zero;

        //planeNormal[partitionAxis] = 1;

        //DebugDraw.DrawPlane(transform.TransformPoint(partitionMeanPoint), planeNormal, extentsMagnitude[partitionAxis]/2f, Color.Lerp(Color.red, Color.blue, (float)depth / 10f), 10 - depth, false);

        // POSITIVE NODE
        Vector3 posMin = parent.bounds.min;
        posMin[partitionAxis] = partitionCoordinate;

        Node posNode         = new Node();
        posNode.bounds       = parentBounds;
        posNode.bounds.min   = posMin;
        parent.positiveChild = posNode;

        
        
        /*
#if DEBUG_KD
        Color c = Color.Lerp(Color.red, Color.blue, (float)depth / 10f);
        c[partitionAxis] = 1f;
        DebugDraw.DrawBounds(parent.positiveChild.bounds, c, 10f - depth, transform);
#endif*/

        // NEGATIVE NODE
        Vector3 negMax = parent.bounds.max;
        negMax[partitionAxis] = partitionCoordinate;

        Node negNode         = new Node();
        negNode.bounds       = parentBounds;
        negNode.bounds.max   = negMax;
        parent.negativeChild = negNode;

        // DebugDraw.DrawBounds(parent.negativeChild.bounds, c, 10f - depth, transform);

        if (positiveTriangles.Count < triCount && positiveTriangles.Count >= 5)
        {
            RecursivePartition(positiveTriangles, depth + 1, posNode, partitionAxis); 
        }
        else
        {
            posNode.triangles = positiveTriangles.ToArray();

            ReturnList(positiveTriangles);

            /*if (drawMeshTreeOnStart)
                DrawTriangleSet(posNode.triangles, DebugDraw.RandomColor(), depth);*/
        }


        if (negativeTriangles.Count < triCount && negativeTriangles.Count >= 5)
        {
            RecursivePartition(negativeTriangles, depth + 1, negNode, partitionAxis);
        }
        else
        {
            negNode.triangles = negativeTriangles.ToArray();

            ReturnList(negativeTriangles);

            /*if (drawMeshTreeOnStart)
                DrawTriangleSet(negNode.triangles, DebugDraw.RandomColor(), depth);*/
        }

    }
    
    void Split(List<int> triangles, float partitionCoordinate, int partitionAxis, List<int> positiveTriangles, List<int> negativeTriangles)
    {
        for (int i=0; i < triangles.Count;i++)
        {
            int triangle = triangles[i];

            bool firstPointAbove  = PointAbovePlane(partitionCoordinate, partitionAxis, vertices[tris[triangle     ]]);
            bool secondPointAbove = PointAbovePlane(partitionCoordinate, partitionAxis, vertices[tris[triangle + 1]]);
            bool thirdPointAbove  = PointAbovePlane(partitionCoordinate, partitionAxis, vertices[tris[triangle + 2]]);

            if (firstPointAbove && secondPointAbove && thirdPointAbove)
            {
                positiveTriangles.Add(triangle);
            }
            else if (!firstPointAbove && !secondPointAbove && !thirdPointAbove)
            {
                negativeTriangles.Add(triangle);
            }
            else
            {
                positiveTriangles.Add(triangle);
                negativeTriangles.Add(triangle);
            }
        }
    }

    //SAH HEURISTICS
    float Heuristic(List<int> triangles, float partitionCoordinate, float partitionStart, float partitionEnd, int partitionAxis, float axisArea, float sideArea)
    {
        int positiveCount = 0;
        int negativeCount = 0;

        for (int i=0; i < triangles.Count; i++)
        {
            int triangle = triangles[i];
            
            bool firstPointAbove  = PointAbovePlane(partitionCoordinate, partitionAxis, vertices[tris[triangle + 0]]);
            bool secondPointAbove = PointAbovePlane(partitionCoordinate, partitionAxis, vertices[tris[triangle + 1]]);
            bool thirdPointAbove  = PointAbovePlane(partitionCoordinate, partitionAxis, vertices[tris[triangle + 2]]);

            if (firstPointAbove && secondPointAbove && thirdPointAbove)
            {
                positiveCount += 1;
            }
            else if (!firstPointAbove && !secondPointAbove && !thirdPointAbove)
            {
                negativeCount += 1;
            }
            else
            {
                //positiveCount += 1;
                //negativeCount += 1;

                //positiveArea += areas[triangle / 3];
                //negativeArea += areas[triangle / 3];
            }
        }

        float ratioLeft = (partitionCoordinate - partitionStart) / (partitionEnd - partitionStart);

        return positiveCount * (sideArea + axisArea * (1f - ratioLeft))
             + negativeCount * (sideArea + axisArea * (ratioLeft));

        /*return positiveCount * Mathf.Sqrt((1f - ratioLeft))
             + negativeCount * Mathf.Sqrt((ratioLeft));*/
    }

    bool PointAbovePlane(Vector3 planeOrigin, Vector3 planeNormal, Vector3 point)
    {
        return Vector3.Dot(point - planeOrigin, planeNormal) >= 0;
    }

    bool PointAbovePlane(float planeCoordinate, int planeAxis, Vector3 point)
    {
        return (point[planeAxis] - planeCoordinate) >= 0;
    }

    float PointDistanceFromPlane(Vector3 planeOrigin, Vector3 planeNormal, Vector3 point)
    {
        return Mathf.Abs(Vector3.Dot((point - planeOrigin), planeNormal));
    }

    // x normal 0 (yz plane)
    // y normal 1 (xz plane)
    // z normal 2 (xy plane)
    float PointDistanceFromPlane(float planeCoordinate, int planeAxis, Vector3 point)
    {
        return Mathf.Abs(point[planeAxis] - planeCoordinate);
    }

    /*float PointDistanceFromPlaneNoAbs(Vector3 planeOrigin, Vector3 planeNormal, Vector3 point)
    {
        return Vector3.Dot(point - planeOrigin, planeNormal);
    }*/

    float PointDistanceFromPlaneNoAbs(float planeCoordinate, int planeAxis, Vector3 point)
    {
        return point[planeAxis] - planeCoordinate;
    }

    /// <summary>
    /// Determines the closest point between a point and a triangle.
    /// Borrowed from RPGMesh class of the RPGController package for Unity, by fholm
    /// The code in this method is copyrighted by the SlimDX Group under the MIT license:
    /// 
    /// Copyright (c) 2007-2010 SlimDX Group
    /// 
    /// Permission is hereby granted, free of charge, to any person obtaining a copy
    /// of this software and associated documentation files (the "Software"), to deal
    /// in the Software without restriction, including without limitation the rights
    /// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    /// copies of the Software, and to permit persons to whom the Software is
    /// furnished to do so, subject to the following conditions:
    /// 
    /// The above copyright notice and this permission notice shall be included in
    /// all copies or substantial portions of the Software.
    /// 
    /// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    /// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    /// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    /// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    /// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    /// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    /// THE SOFTWARE.
    /// 
    /// </summary>
    /// <param name="point">The point to test.</param>
    /// <param name="vertex1">The first vertex to test.</param>
    /// <param name="vertex2">The second vertex to test.</param>
    /// <param name="vertex3">The third vertex to test.</param>
    /// <param name="result">When the method completes, contains the closest point between the two objects.</param>
    void ClosestPointOnTriangleToPoint(ref Vector3 vertex1, ref Vector3 vertex2, ref Vector3 vertex3, ref Vector3 point, out Vector3 result)
    {
        //Source: Real-Time Collision Detection by Christer Ericson
        //Reference: Page 136

        //Check if P in vertex region outside A
        Vector3 ab = vertex2 - vertex1;
        Vector3 ac = vertex3 - vertex1;
        Vector3 ap = point - vertex1;

        float d1 = Vector3.Dot(ab, ap);
        float d2 = Vector3.Dot(ac, ap);
        if (d1 <= 0.0f && d2 <= 0.0f)
        {
            result = vertex1; //Barycentric coordinates (1,0,0)
            return;
        }

        //Check if P in vertex region outside B
        Vector3 bp = point - vertex2;
        float d3 = Vector3.Dot(ab, bp);
        float d4 = Vector3.Dot(ac, bp);
        if (d3 >= 0.0f && d4 <= d3)
        {
            result = vertex2; // barycentric coordinates (0,1,0)
            return;
        }

        //Check if P in edge region of AB, if so return projection of P onto AB
        float vc = d1 * d4 - d3 * d2;
        if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f)
        {
            float v = d1 / (d1 - d3);
            result = vertex1 + v * ab; //Barycentric coordinates (1-v,v,0)
            return;
        }

        //Check if P in vertex region outside C
        Vector3 cp = point - vertex3;
        float d5 = Vector3.Dot(ab, cp);
        float d6 = Vector3.Dot(ac, cp);
        if (d6 >= 0.0f && d5 <= d6)
        {
            result = vertex3; //Barycentric coordinates (0,0,1)
            return;
        }

        //Check if P in edge region of AC, if so return projection of P onto AC
        float vb = d5 * d2 - d1 * d6;
        if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f)
        {
            float w = d2 / (d2 - d6);
            result = vertex1 + w * ac; //Barycentric coordinates (1-w,0,w)
            return;
        }

        //Check if P in edge region of BC, if so return projection of P onto BC
        float va = d3 * d6 - d5 * d4;
        if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f)
        {
            float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            result = vertex2 + w * (vertex3 - vertex2); //Barycentric coordinates (0,1-w,w)
            return;
        }

        //P inside face region. Compute Q through its barycentric coordinates (u,v,w)
        float denom = 1.0f / (va + vb + vc);
        float v2 = vb * denom;
        float w2 = vc * denom;
        result = vertex1 + ab * v2 + ac * w2; //= u*vertex1 + v*vertex2 + w*vertex3, u = va * denom = 1.0f - v - w
    }

    void DrawTriangleSet(int[] triangles, Color color, float duration, Transform transform)
    {
        foreach (int triangle in triangles)
        {
            DebugDraw.DrawTriangle(vertices[tris[triangle]], vertices[tris[triangle + 1]], vertices[tris[triangle + 2]], color, transform, duration);
        }
    }

    public void RecursiveDraw(Node node, int depth, bool permutation)
    {
        var color = Color.HSVToRGB(((permutation?0.05f:0f) + depth * 0.05f) % 1f, 1f, 1f);
        color.a = Mathf.Clamp((depth)/10f, 0f, 1f) * 0.1f;
        Gizmos.color = color;

        Gizmos.DrawWireCube(node.bounds.center, node.bounds.size);

        if (node.positiveChild != null)
        {
            RecursiveDraw(node.positiveChild, depth + 1, permutation);
        }

        if (node.negativeChild != null)
        {
            RecursiveDraw(node.negativeChild, depth + 1, !permutation);
        }
    }

}