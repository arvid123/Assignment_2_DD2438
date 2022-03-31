using Priority_Queue;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;
using System.Diagnostics;
using System.IO;
using Priority_Queue;

namespace Assets.Scrips
{
    class Pathgen
    {
        TerrainManager terrain_manager;
        float terrain_padding;
        List<Vector3> my_path;
        Stack<Waypoint> optimal_path;
        float max_turning_velocity;
        public Vector3[] wps;
        List<Vector3> friend_list = new List<Vector3>();
        List<Vector3> enemy_list = new List<Vector3>();
        bool[,] neighbors;
        int nr_wps;
        int nr_friends;
        int nr_enemies;


        public Pathgen(TerrainManager t, float tp, float mtv, string vehicle, List<Vector3> enemies, List<Vector3> friends)
        {
            terrain_manager = t;
            terrain_padding = tp;
            max_turning_velocity = mtv;

            //Vector3 start_pos = terrain_manager.myInfo.start_pos;
            //Vector3 goal_pos = terrain_manager.myInfo.goal_pos;

            my_path = new List<Vector3>();
            //UnityEngine.Debug.DrawRay(start_pos, Vector3.forward, Color.magenta, 1000f);

            friend_list = friends;
            enemy_list = enemies;


            // Plan your path here
            // ...

            //for (int i = 0; i < 3; i++)
            //{
            //}


            // Create padded visibility graph
            List<GameObject> padded_cubes = new List<GameObject>();
            var ter = terrain_manager.myInfo;
            float x_step = (ter.x_high - ter.x_low) / ter.x_N;
            float z_step = (ter.z_high - ter.z_low) / ter.z_N;
            for (int i = 0; i < ter.x_N; i++)
            {
                for (int j = 0; j < ter.z_N; j++)
                {
                    if (ter.traversability[i, j] > 0.5f)
                    {
                        // Add invisible padded cube for collision detection
                        GameObject cube = GameObject.CreatePrimitive(PrimitiveType.Cube);
                        cube.transform.position = new Vector3(ter.get_x_pos(i), 0.0f, ter.get_z_pos(j));
                        cube.transform.localScale = new Vector3(x_step + terrain_padding * 1.9f, 15.0f, z_step + terrain_padding * 1.9f);
                        cube.GetComponent<MeshRenderer>().enabled = false;
                        cube.layer = 1;
                        padded_cubes.Add(cube);

                        Vector3 lower_left = new Vector3(ter.get_x_pos(i) - x_step / 2 - terrain_padding, 0, ter.get_z_pos(j) - z_step / 2 - terrain_padding);
                        if (ter.traversability[ter.get_i_index(lower_left.x), ter.get_j_index(lower_left.z)] <= 0.5f)
                        {
                            // Make sure we only add convex corners to the graph
                            if (j > 0 && ter.traversability[i, j - 1] <= 0.5f && i > 0 && ter.traversability[i - 1, j] <= 0.5f)
                            {
                                my_path.Add(lower_left);
                            }
                        }
                        Vector3 lower_right = new Vector3(ter.get_x_pos(i) + x_step / 2 + terrain_padding, 0, ter.get_z_pos(j) - z_step / 2 - terrain_padding);
                        if (ter.traversability[ter.get_i_index(lower_right.x), ter.get_j_index(lower_right.z)] <= 0.5f)
                        {
                            // Make sure we only add convex corners to the graph
                            if (j > 0 && ter.traversability[i, j - 1] <= 0.5f && i < ter.x_N - 1 && ter.traversability[i + 1, j] <= 0.5f)
                            {
                                my_path.Add(lower_right);
                            }
                        }
                        Vector3 upper_left = new Vector3(ter.get_x_pos(i) - x_step / 2 - terrain_padding, 0, ter.get_z_pos(j) + z_step / 2 + terrain_padding);
                        if (ter.traversability[ter.get_i_index(upper_left.x), ter.get_j_index(upper_left.z)] <= 0.5f)
                        {
                            // Make sure we only add convex corners to the graph
                            if (j < ter.z_N - 1 && ter.traversability[i, j + 1] <= 0.5f && i > 0 && ter.traversability[i - 1, j] <= 0.5f)
                            {
                                my_path.Add(upper_left);
                            }
                        }
                        Vector3 upper_right = new Vector3(ter.get_x_pos(i) + x_step / 2 + terrain_padding, 0, ter.get_z_pos(j) + z_step / 2 + terrain_padding);
                        if (ter.traversability[ter.get_i_index(upper_right.x), ter.get_j_index(upper_right.z)] <= 0.5f)
                        {
                            // Make sure we only add convex corners to the graph
                            if (j < ter.z_N - 1 && ter.traversability[i, j + 1] <= 0.5f && i < ter.x_N - 1 && ter.traversability[i + 1, j] <= 0.5f)
                            {
                                my_path.Add(upper_right);
                            }
                        }
                    }
                }
            }
            nr_wps = friend_list.Count + enemy_list.Count + my_path.Count;
            wps = new Vector3[nr_wps];
            nr_friends = friend_list.Count;
            nr_enemies = enemy_list.Count;

            // Add all friends to wps
            for (int i = 0; i < friend_list.Count; i++)
            {
                wps[i] = friend_list[i];
            }
            // Add all enemies to wps
            for (int i = 0; i < enemy_list.Count; i++)
            {
                wps[friend_list.Count + i] = enemy_list[i];
            }
            // Add my_path (visibility graph) to wps
            for (int i = 0; i < my_path.Count; i++)
            {
                wps[friend_list.Count + enemy_list.Count + i] = my_path[i];
            }

            // neighbors[x, y] == true means x, y are neighbors and also neighbors[y, x] == true
            neighbors = new bool[nr_wps, nr_wps];
            
            for (int i = 0; i < nr_wps; i++)
            {
                for (int j = 0; j < nr_wps; j++)
                {
                    if (i != j && !Physics.Linecast(wps[i], wps[j], LayerMask.GetMask("TransparentFX")))
                    {
                        neighbors[i, j] = true;
                        neighbors[j, i] = true;
                        //UnityEngine.Debug.DrawLine(w.pos, otherw.pos, Color.red, 100f);

               
                    }
                    // i is an enemy
                    else if (i != j && (i < nr_friends + nr_enemies && i >= nr_friends))
                    {
                        bool enemy_within_padding = Physics.OverlapSphere(wps[i], 0, LayerMask.GetMask("TransparentFX")).Length > 0;

                        if(enemy_within_padding && !Physics.Linecast(wps[i], wps[j], LayerMask.GetMask("CubeWalls"))){
                            neighbors[i, j] = true;
                            neighbors[j, i] = true;
                            //UnityEngine.Debug.DrawLine(w.pos, otherw.pos, Color.red, 100f
                                
                        }
                    }
                }
            }

            // Delete padded cubes
            foreach (var cube in padded_cubes)
            {
                cube.SetActive(false);
            }
        }

        public Stack<Waypoint> getOptimalPath()
        {
            return new Stack<Waypoint>(new Stack<Waypoint>(optimal_path));
        }

        public List<Vector3> getOptimalPathList()
        {
            return new List<Waypoint>(optimal_path).ConvertAll(x => x.pos);
        }

        public List<List<Vector3>> vrpSolver()
        {

            // Run vrp solver
            // Use ProcessStartInfo class
            ProcessStartInfo startInfo = new ProcessStartInfo();
            startInfo.WorkingDirectory = Directory.GetCurrentDirectory() + "\\Assets\\Resources\\or-tools\\bin";
            startInfo.CreateNoWindow = false;
            startInfo.UseShellExecute = true;
            startInfo.FileName = "vrp_new.exe";
            startInfo.WindowStyle = ProcessWindowStyle.Hidden;
            UnityEngine.Debug.Log(Directory.GetCurrentDirectory());

            startInfo.Arguments += String.Format("{0}", nr_enemies + 1);
            // Create distance matrix
            for (int i = nr_friends - 1; i < nr_friends + nr_enemies; i++)
            {
                for (int j = nr_friends - 1; j < nr_friends + nr_enemies; j++)
                {
                    startInfo.Arguments += " ";
                    startInfo.Arguments += pathCost(i, j);
                    
                }
            }

            try
            {
                // Start the process with the info we specified.
                // Call WaitForExit and then the using statement will close.
                using (Process exeProcess = Process.Start(startInfo))
                {
                    exeProcess.WaitForExit();
                }
            }
            catch
            {
                // Log error.
            }

            List<int>[] vrp_routes = new List<int>[nr_friends];

            string[] vrp_solutions = System.IO.File.ReadAllLines(Directory.GetCurrentDirectory() + "\\Assets\\Resources\\or-tools\\bin\\solution.txt");
            for (int i = 0; i < nr_friends; i++)
            {
                List<int> route = vrp_solutions[i].Split(' ').Select(Int32.Parse).ToList();
                route.RemoveAt(route.Count - 1);
                route[0] = -i;
                UnityEngine.Debug.Log(string.Join(" ", route.ToArray()));
                vrp_routes[i] = route;
            }

            List<List<Vector3>> smooth_paths = new List<List<Vector3>>();
            for (int i = 0; i < nr_friends; i++)
            {
                List<Vector3> smooth_path = new List<Vector3>();
                for (int j = 1; j < vrp_routes[i].Count; j++)
                {
                    smooth_path.AddRange(getBezierPathList(vrp_routes[i][j - 1] + nr_friends-1, vrp_routes[i][j] + nr_friends-1));
                }
                smooth_paths.Add(smooth_path);
            }

            return smooth_paths;
        }

        public List<Vector3> tspSolver(int[] targets)
        {
            List<int> path = new List<int>();
            var start = targets[0];
            Dictionary<int, bool> marked = new Dictionary<int, bool>();
            for (int i = 0; i < targets.Length; i++)
            {
                marked.Add(targets[i], false);
            }
            path.Add(start);
            marked[start] = true;

            while (true)
            {
                int best_target = start;
                float lowest_cost = float.PositiveInfinity;

                foreach (var target in targets)
                {
                    if (!marked[target])
                    {
                        float cost = pathCost(start, target);
                        if (cost < lowest_cost)
                        {
                            lowest_cost = cost;
                            best_target = target;
                            start = target;
                        }
                    }
                }
                marked[best_target] = true;
                path.Add(best_target);

                // All targets have been marked
                if (!marked.ContainsValue(false))
                {
                    break;
                }
            }

            List<Vector3> smooth_path = new List<Vector3>();
            for(int i = 1; i < path.Count; i++)
            {
                smooth_path.AddRange(getBezierPathList(path[i-1], path[i]));
            }

            return smooth_path;
        }

        public float pathCost(int start, int goal)
        {
            var path = A_star(start, goal, car_g, car_h);
            float sum_g = 0;
            Vector3 old_wp = path[0];
            foreach(var wp in path)
            {
                sum_g += Vector3.Distance(old_wp, wp);
                old_wp = wp;
            }
            return sum_g;
        }

        // https://en.wikipedia.org/wiki/A*_search_algorithm 
        public List<Vector3> A_star(int start, int goal, Func<int, int, float> g, Func<int, int, float> h)
        {
            SimplePriorityQueue<int, double> openSet = new SimplePriorityQueue<int, double>();
            openSet.Enqueue(start, h(start, goal));

            int[] cameFrom = new int[nr_wps];

            // g score map
            float[] g_scores = new float[nr_wps];
            for (int i = 0; i < g_scores.Length; i++)
            {
                g_scores[i] = float.PositiveInfinity;
            }
            g_scores[start] = 0;

            while (openSet.Count > 0)
            {
                int current = openSet.Dequeue();
                if (current == goal)
                {
                    return reconstruct_path(cameFrom, start, goal);
                }

                foreach (int neighbor in getNeighborsList(current))
                {
                    //UnityEngine.Debug.DrawLine(current.pos, neighbor.pos, Color.green, 100f);
                    float tentative_gScore = g_scores[current] + g(current, neighbor);

                    if (tentative_gScore < g_scores[neighbor])
                    {
                        // This path to neighbor is better than any other recorded one.
                        cameFrom[neighbor] = current;
                        g_scores[neighbor] = tentative_gScore;
                        openSet.EnqueueWithoutDuplicates(neighbor, tentative_gScore + h(neighbor, goal));
                    }
                }
            }

            // Should be unreachable
            UnityEngine.Debug.Log("A* failed");
            return null;
        }

        private List<int> getNeighborsList(int wp)
        {
            var neighbors_list = new List<int>();
            for (int i = 0; i < nr_wps; i++)
            {
                if (neighbors[i, wp])
                {
                    neighbors_list.Add(i);
                }
            }
            return neighbors_list;
        }

        private float car_g(int current, int neighbor)
        {
            float g = Vector3.Distance(wps[current], wps[neighbor]);
            return g;
        }

        private double drone_g(Waypoint start, Waypoint goal, Waypoint current, Waypoint neighbor, List<Waypoint> wps)
        {
            double g = Vector3.Distance(current.pos, neighbor.pos);

            if (current.cameFrom != null)
            {
                g += Vector3.Angle(current.cameFrom.pos - current.pos, neighbor.pos - current.pos) * 0f;
            }

            return g;
        }

        // The heuristic for the A* algorithm
        private double drone_h(Waypoint w, Waypoint goal)
        {
            return Vector3.Distance(w.pos, goal.pos);
        }

        private float car_h(int w, int goal)
        {
            return Vector3.Distance(wps[w], wps[goal]);
        }

        // Create stack of waypoints and also add drone_goal_vel depending on the angle between the points
        private List<Vector3> reconstruct_path(int[] cameFrom, int start, int goal)
        {
            List<Vector3> path = new List<Vector3>();
            int current = goal;
            while (current != start)
            {
                path.Add(wps[current]);
                current = cameFrom[current];
            }
            path.Add(wps[start]);
            path.Reverse();
            return path;

        }

        private float squared(float x)
        {
            return x * x;
        }

        private float cubed(float x)
        {
            return x * x * x;
        }


        public Stack<Waypoint> getBezierPath()
        {
            var coarse_path = new List<Waypoint>(getOptimalPath());
            var bezier_path = new List<Waypoint>();

            float slide = 6f;

            bezier_path.Add(coarse_path[0]);
            for (int i = 1; i < coarse_path.Count - 1; i++)
            {
                Vector3 backwards_direction = (coarse_path[i - 1].pos - coarse_path[i].pos).normalized;
                Vector3 forwards_direction = (coarse_path[i + 1].pos - coarse_path[i].pos).normalized;

                Vector3 control_point_1 = coarse_path[i].pos + backwards_direction * Mathf.Min(slide, (coarse_path[i - 1].pos - coarse_path[i].pos).magnitude / 2);
                Vector3 control_point_3 = coarse_path[i].pos + forwards_direction * Mathf.Min(slide, (coarse_path[i + 1].pos - coarse_path[i].pos).magnitude / 2);

                for (float t = 0; t <= 1; t += 0.02f)
                {
                    Vector3 sp = BezierCurve.Quadratic(new Vector2(control_point_1.x, control_point_1.z), new Vector2(coarse_path[i].pos.x, coarse_path[i].pos.z), new Vector2(control_point_3.x, control_point_3.z), t);
                    Vector3 sample_point = new Vector3(sp.x, 0, sp.y);
                    bezier_path.Add(new Waypoint(sample_point));
                }
            }
            bezier_path.Add(coarse_path[coarse_path.Count - 1]);

            for (int i = 1; i < bezier_path.Count - 1; i++)
            {
                float turn_angle_ratio = Mathf.Pow(Vector3.Angle(bezier_path[i - 1].pos - bezier_path[i].pos, bezier_path[i + 1].pos - bezier_path[i].pos) / 180f, 12);
                bezier_path[i].drone_goal_vel = (bezier_path[i + 1].pos - bezier_path[i].pos).normalized * max_turning_velocity * turn_angle_ratio;
            }
            for (int i = 0; i < bezier_path.Count - 25; i++)
            {
                if (Vector3.Distance(bezier_path[i + 25].pos, bezier_path[i].pos) < 10f && bezier_path[i + 25].drone_goal_vel.magnitude < bezier_path[i].drone_goal_vel.magnitude)
                {
                    bezier_path[i].drone_goal_vel = bezier_path[i].drone_goal_vel.normalized * Mathf.Lerp(bezier_path[i].drone_goal_vel.magnitude, bezier_path[i + 25].drone_goal_vel.magnitude, 0.8f);
                }
            }

            bezier_path[0].drone_goal_vel = (bezier_path[1].pos - bezier_path[0].pos).normalized * 15f;
            bezier_path[bezier_path.Count - 1].drone_goal_vel = (bezier_path[bezier_path.Count - 1].pos - bezier_path[bezier_path.Count - 2].pos).normalized * 15f;

            // draw goal_vels
            foreach (Waypoint w in bezier_path)
            {
                ////UnityEngine.Debug.DrawLine(w.pos, w.pos + w.drone_goal_vel, Color.yellow, 1000f);
            }

            return new Stack<Waypoint>(new Stack<Waypoint>(bezier_path));
        }

        public List<Vector3> getBezierPathList(int start, int goal)
        {
            var coarse_path = A_star(start, goal, car_g, car_h);
            var bezier_path = new List<Vector3>();

            float slide = 6f;

            bezier_path.Add(coarse_path[0]);
            for (int i = 1; i < coarse_path.Count - 1; i++)
            {
                Vector3 backwards_direction = (coarse_path[i - 1] - coarse_path[i]).normalized;
                Vector3 forwards_direction = (coarse_path[i + 1] - coarse_path[i]).normalized;

                Vector3 control_point_1 = coarse_path[i] + backwards_direction * Mathf.Min(slide, (coarse_path[i - 1] - coarse_path[i]).magnitude / 2);
                Vector3 control_point_3 = coarse_path[i] + forwards_direction * Mathf.Min(slide, (coarse_path[i + 1] - coarse_path[i]).magnitude / 2);

                for (float t = 0; t <= 1; t += 0.02f)
                {
                    Vector3 sp = BezierCurve.Quadratic(new Vector2(control_point_1.x, control_point_1.z), new Vector2(coarse_path[i].x, coarse_path[i].z), new Vector2(control_point_3.x, control_point_3.z), t);
                    Vector3 sample_point = new Vector3(sp.x, 0, sp.y);
                    bezier_path.Add(sample_point);
                }
            }
            bezier_path.Add(coarse_path[coarse_path.Count - 1]);
            return bezier_path;
        }


        public Stack<Waypoint> getSmoothPath()
        {
            var chosen_path = getOptimalPath();
            var chosen_path_list = new List<Waypoint>(new Stack<Waypoint>(getOptimalPath()));
            var coarse_path = new List<Waypoint>(chosen_path);
            List<Vector3> smooth_path;
            // linear interpolation
            int chosen_path_len = 1;
            int cps_len = 1;
            int num_interpolation;
            Waypoint current = chosen_path.Pop();
            List<Vector3> cps = new List<Vector3>(); // create contol points
            cps.Add(current.pos);
            List<Vector3> cps_new = new List<Vector3>();
            int cps_new_len = 1;
            cps_new.Add(current.pos);
            // UnityEngine.UnityEngine.Debug.Log("cps" + current.pos);

            while (chosen_path.Count > 1)
            {
                Waypoint next = chosen_path.Pop();
                float linear_dis = Vector3.Distance(current.pos, next.pos);
                //UnityEngine.UnityEngine.Debug.Log("Linear distance" + linear_dis);
                //num_interpolation = (int)Math.Ceiling(linear_dis / 2);
                num_interpolation = 11; // 
                //UnityEngine.UnityEngine.Debug.Log("num_interpolation" + num_interpolation);
                for (int i = 1; i <= num_interpolation; i++)
                {
                    float rate = (float)i / num_interpolation;
                    Vector3 add_point = Vector3.Lerp(current.pos, next.pos, rate);
                    cps.Add(add_point);
                    //UnityEngine.UnityEngine.Debug.Log("cps"+ add_point);
                    if (i == 2) // //
                        cps_new.Add(add_point); // //
                    if (i == 9) // //
                        cps_new.Add(add_point); // //
                }
                cps_new_len += 2;
                // cps_new.Add(next.pos);// //
                // cps_new_len++; // //
                chosen_path_len++;
                cps_len += num_interpolation;
                current = next;
            }
            cps_new.Add(current.pos);
            cps_new_len += 1;
            //UnityEngine.UnityEngine.Debug.Log("Linear Path Length " + chosen_path_len);
            //UnityEngine.UnityEngine.Debug.Log("Possible Control points Length " + cps_len);

            /*
            List<Vector3> cps_new = new List<Vector3>(); // create contol points
            int cps_new_len = 0;
            int interval = 4;
            for (int i = 0; i < cps_len - 1; i++)
            {
                if ((i % interval) == 0)
                {
                    cps_new.Add(cps[i]);
                    //UnityEngine.UnityEngine.Debug.Log("cps_new" + cps[i]);
                    cps_new_len++;
                }
            }
            cps_new.Add(cps[cps_len - 1]);
            //UnityEngine.UnityEngine.Debug.Log("cps_new" + cps[cps_len - 1]);
            cps_new_len++;
            //UnityEngine.UnityEngine.Debug.Log("Control points Length" + cps_new_len);
            
            */
            Vector3 old_cp = cps_new[0];
            foreach (var cp in cps_new)
            {
                ////UnityEngine.Debug.DrawLine(old_cp, cp, Color.red, 100f);
                old_cp = cp;
            }


            // Smoothness
            SplineCurve curve = new SplineCurve();
            foreach (var cp in cps_new)
            {
                curve.AddNode(cp);
            }
            curve.AddCatmull_RomControl();
            smooth_path = new List<Vector3>(); //create smooth path
            for (int i = 0; i < curve.segmentList.Count; i++)
            {
                float add = 1f / 10;  // 表示两个关键点之间取20个点，可根据需要设置
                for (float j = 0; j < 1; j += add)
                {
                    Vector3 point = curve.segmentList[i].GetPoint(j);
                    smooth_path.Add(point);
                }
            }

            Vector3 old_sp = smooth_path[0];
            foreach (var sp in smooth_path)
            {
                UnityEngine.Debug.DrawLine(old_sp, sp, Color.blue, 100f);
                old_sp = sp;
            }
            //UnityEngine.UnityEngine.Debug.Log("Smooth Path Length " + smooth_path_len);

            Stack<Waypoint> path = new Stack<Waypoint>(smooth_path.ConvertAll<Waypoint>(x => new Waypoint(x)));
            Waypoint current_point = path.Pop();
            Waypoint next_point = path.Pop();
            Waypoint next_next = path.Pop();
            Stack<Waypoint> smooth_stack = new Stack<Waypoint>();

            do
            {

                float turn_angle_ratio = Mathf.Pow(Vector3.Angle(current_point.pos - next_point.pos, next_next.pos - next_point.pos) / 180f, 10);
                next_point.drone_goal_vel = -(next_next.pos - next_point.pos).normalized * max_turning_velocity * turn_angle_ratio;

                // TODO: Assign drone_goal_vel of closest point in coarse_path

                //next_point.drone_goal_vel = -(next_next.pos - next_point.pos).normalized * closestPoint(next_point, coarse_path).drone_goal_vel.magnitude;

                smooth_stack.Push(current_point);
                current_point = next_point;
                next_point = next_next;
                next_next = path.Pop();
            } while (path.Count > 0);
            smooth_stack.Push(next_next);

            return smooth_stack;
        }

        private Waypoint closestPoint(Waypoint point, List<Waypoint> coarse_path)
        {
            Waypoint closest = coarse_path[0];

            foreach (Waypoint wp in coarse_path)
            {
                if (Vector3.Distance(wp.pos, point.pos) < Vector3.Distance(closest.pos, point.pos))
                {
                    closest = wp;
                }
            }

            return closest;
        }
    }
}
