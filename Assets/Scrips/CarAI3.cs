using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using Assets.Scrips;

namespace UnityStandardAssets.Vehicles.Car
{
    [RequireComponent(typeof(CarController))]
    public class CarAI3 : MonoBehaviour
    {
        private CarController m_Car; // the car controller we want to use

        public GameObject terrain_manager_game_object;
        TerrainManager terrain_manager;

        public GameObject[] friends;
        public GameObject[] enemies;
        List<Vector3> route { get; set; }

        private void Start()
        {
            // get the car controller
            m_Car = GetComponent<CarController>();
            terrain_manager = terrain_manager_game_object.GetComponent<TerrainManager>();


            // note that both arrays will have holes when objects are destroyed
            // but for initial planning they should work
            friends = GameObject.FindGameObjectsWithTag("Player");
            enemies = GameObject.FindGameObjectsWithTag("Enemy");
            foreach (GameObject obj in enemies)
            {
                Debug.DrawLine(transform.position, obj.transform.position, Color.black, 10f);
            }

            route = new List<Vector3>();

            // Plan your path here
            // ...
            int i = 0;
            route.Add(transform.position);

            foreach(var enemy in enemies)
            {
                friends[i].GetComponent<CarAI3>().route.Add(enemy.transform.position);
                i = (i + 1) % 3;
            }
            Debug.Log(friends[0].GetComponent<CarAI3>().route[0]);
            Debug.Log(friends[1].GetComponent<CarAI3>().route[0]);
            Debug.Log(friends[2].GetComponent<CarAI3>().route[0]);

            Pathgen pg = new Pathgen(terrain_manager, 4f, 5f, "car", enemies, friends);
            List<Vector3> good_route = new List<Vector3>();

            for(i = 0; i < route.Count-1; i++)
            {
                pg = new Pathgen(terrain_manager, 4f, 5f, "car", enemies, friends);
                var segment = pg.getBezierPathList(route[i], route[i + 1]);
                good_route.AddRange(segment);
            }
            
            //Debug.Log(String.Format("Goal at {0}, {1}", route[1].x, route[1].z));
            
            var old = good_route[0];
            Debug.Log(good_route.Count);
            foreach(Vector3 pos in good_route)
            {
                Debug.Log("Long");
                Debug.DrawLine(old, pos, Color.blue, 1000f);
                old = pos;
            }
        }


        private void FixedUpdate()
        {


            // Execute your path here
            // ...

            Vector3 avg_pos = Vector3.zero;

            foreach (GameObject friend in friends)
            {
                avg_pos += friend.transform.position;
            }
            avg_pos = avg_pos / friends.Length;
            Vector3 direction = (avg_pos - transform.position).normalized;

            bool is_to_the_right = Vector3.Dot(direction, transform.right) > 0f;
            bool is_to_the_front = Vector3.Dot(direction, transform.forward) > 0f;

            float steering = 0f;
            float acceleration = 0;

            if (is_to_the_right && is_to_the_front)
            {
                steering = 1f;
                acceleration = 1f;
            }
            else if (is_to_the_right && !is_to_the_front)
            {
                steering = -1f;
                acceleration = -1f;
            }
            else if (!is_to_the_right && is_to_the_front)
            {
                steering = -1f;
                acceleration = 1f;
            }
            else if (!is_to_the_right && !is_to_the_front)
            {
                steering = 1f;
                acceleration = -1f;
            }

            // this is how you access information about the terrain
            int i = terrain_manager.myInfo.get_i_index(transform.position.x);
            int j = terrain_manager.myInfo.get_j_index(transform.position.z);
            float grid_center_x = terrain_manager.myInfo.get_x_pos(i);
            float grid_center_z = terrain_manager.myInfo.get_z_pos(j);

            Debug.DrawLine(transform.position, new Vector3(grid_center_x, 0f, grid_center_z));


            // this is how you control the car
            Debug.Log("Steering:" + steering + " Acceleration:" + acceleration);
            m_Car.Move(steering, acceleration, acceleration, 0f);
            //m_Car.Move(0f, -1f, 1f, 0f);


        }
    }
}
