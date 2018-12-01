#include "hlt/game.hpp"
#include "hlt/constants.hpp"
#include "hlt/log.hpp"

#include <random>
#include <ctime>
#include <fstream>
#include <math.h>

using namespace std;
using namespace hlt;

struct Node
{
	int halite;
	Position position;
	int id;
	int fleet = -1;
};

struct FleetShip
{
	shared_ptr<Ship> ship;
	Position destination;
	int row;
	int state;
	int fleet;
};

struct Fleet
{
	FleetShip ships[4];
	int id;
	vector<int> nodes;
};

vector<Node> nodes;
vector<Fleet> fleets;
int minHal = 0;

int nextNode(FleetShip &ship, int DIM)
{
    Position pos = ship.ship->pos;
    int currentNode = (pos.x / 4) + (pos.y / 4) * (DIM / 4);
    for (int i : fleets[ship.fleet].nodes)
    {
        if (i == currentNode)
            return i;
    }
    int[] nextNodes = {currentNode + 1, currentNode - 1, currentNode - DIM/4, currentNode + DIM/4};
    sort(nextNodes.begin(), nextNodes.end, bool func(int a, int b) { return nodes[a].halite > nodes[b].halite});
    for (int i : nextNodes)
    {
        if (nodes[i].fleet == -1)
        {
            fleets[ship.fleet].push_back(i);
            nodes[i].fleet = ship.fleet;
            return i;
        }
    }
    return -1;
}

void addToFleet(shared_ptr<Ship> ship)
{
    if (fleets.size() == 0 || fleets.back().ships.size() == 4)
    {
        for (unsigned int i = 0; i < nodes.size(); ++i)
        {
            int initNode = -1;
            if (nodes[nextNodes].fleet == -1)
            {
                initNode = i;
                nodes[i].fleet = fleets.size();
                break;
            }
        }
        fleets.push_back(Fleet{{}, fleets.size(), {initNode}});
    }
	fleets.back().ships.push_back(FleetShip{ship, Position(nodes[i].position.x, nodes[i].y + fleets.back().ships.size(), fleets.back().ships.size(), 1, fleets.back().id});
}

Position updateDestination(FleetShip ship, int DIM, shared_ptr<Player> me) {
	int posx = ship.ship->pos.x;
	int posy = ship.ship->pos.y;
	int currentNode = (pos.x / 4) + (pos.y / 4) * (DIM / 4);

	if (x % 4 == 3 && ship.state == 2) {
		Node n = nodes[nextNode(ship, DIM)];
		ship.destination = Position(n.position.x, n.position.y + posy);
		ship.state == 3;
	}
	else if (ship.state == 2){
		mine(ship);
	}
	else if (ship.ship->halite == 1000) {
		ship.state == 0;
		ship.destination = me->shipyard->position;
	}
	else if (ship.state == 3 && ship.destination == ship.ship->pos) {
		ship.state == 2;
	}
}

void mine(FleetShip ship, Game game) {
	if (game.game_map.at(Position(i, j))->halite < minHal) {
		command_queue.push_back(ship.ship->move(Direction.East));
	}
	else {
		command_queue.push_back(ship.ship->stay_still());
	}
}

void moveShip(FleetShip ship) {
	if (x % 4 != 3 && ship.state != 2) {
		//move ship
	}
}

/*
int main(int argc, char* argv[]) {
	unsigned int rng_seed;
	if (argc > 1) {
		rng_seed = static_cast<unsigned int>(stoul(argv[1]));
	}
	else {
		rng_seed = static_cast<unsigned int>(time(nullptr));
	}
	mt19937 rng(rng_seed);

	Game game;

	int DIM = game.game_map.cells.size();
	for (unsigned int i = 0; i < DIM; i += 4)
	{
		for (unsigned int j = 0; j < DIM; j += 4)
		{
			nodes.push_back(Node{ 0, Position(i, j), (i / 4 + j / 4 * DIM / 4) };
		}
	}
	int totHal = 0;
	for (unsigned int i = 0; i < DIM; ++i)
	{
		for (unsigned int j = 0; j < DIM; ++j)
		{
			totHal += game.game_map.at(Position(i, j))->halite;
			nodes[i / 4 + j / 4 * DIM / 4].halite += game.game_map.at(Position(i, j))->halite;
		}
	}

	int meanHal = totHal / (DIM*DIM);
	int stdDev = 0;

	for (unsigned int i = 0; i < DIM; ++i)
	{
		for (unsigned int j = 0; j < DIM; ++j)
		{
			stdDev += (game.game_map.at(Position(i, j))->halite - meanHal)(game.game_map.at(Position(i, j))->halite - meanHal);
		}
	}

	stdDev /= (DIM*DIM);
	stdDev = sqrt(stDEV);
	minHal = meanHal - 2 * stdDev;

	// At this point "game" variable is populated with initial map data.
	// This is a good place to do computationally expensive start-up pre-processing.
	// As soon as you call "ready" function below, the 2 second per turn timer will start.

	game.ready("MyCppBot");

	log::log("Successfully created bot! My Player ID is " + to_string(game.my_id) + ". Bot rng seed is " + to_string(rng_seed) + ".");
	int init = 0;
	for (;;) {
		game.update_frame();
		shared_ptr<Player> me = game.me;
		unique_ptr<GameMap>& game_map = game.game_map;

		vector<Command> command_queue;

		for (const auto& ship_iterator : me->ships) {
			shared_ptr<Ship> ship = ship_iterator.second;

			bool added = false;
            		for (EntityId id : addedShips)
            		{
                		if (id == ship.id)
                    		added = true;
            		}
            		if (!added)
            		{
                		addToFleet(ship);
            		}
			
			if (game_map->at(ship)->halite < constants::MAX_HALITE / 10 || ship->is_full()) {
				Direction random_direction = ALL_CARDINALS[rng() % 4];
				command_queue.push_back(ship->move(random_direction));
			}
			else {
				command_queue.push_back(ship->stay_still());
			}
		}

		if (
			game.turn_number <= 200 &&
			me->ships.size() < 4 &&
			me->halite >= constants::SHIP_COST &&
			!game_map->at(me->shipyard)->is_occupied())
		{
			command_queue.push_back(me->shipyard->spawn());
		}

		if (!game.end_turn(command_queue)) {
			break;
		}
	}

	return 0;
}
*/

int main(int argc, char* argv[]) {
	unsigned int rng_seed;
	if (argc > 1) {
		rng_seed = static_cast<unsigned int>(stoul(argv[1]));
	}
	else {
		rng_seed = static_cast<unsigned int>(time(nullptr));
	}
	mt19937 rng(rng_seed);

	Game game;
	// At this point "game" variable is populated with initial map data.
	// This is a good place to do computationally expensive start-up pre-processing.
	// As soon as you call "ready" function below, the 2 second per turn timer will start.
	game.ready("MyCppBot");

	log::log("Successfully created bot! My Player ID is " + to_string(game.my_id) + ". Bot rng seed is " + to_string(rng_seed) + ".");

	for (;;) {
		game.update_frame();
		shared_ptr<Player> me = game.me;
		unique_ptr<GameMap>& game_map = game.game_map;

		vector<Command> command_queue;

		for (const auto& ship_iterator : me->ships) {
			shared_ptr<Ship> ship = ship_iterator.second;
			if (game_map->at(ship)->halite < constants::MAX_HALITE / 10 || ship->is_full()) {
				Direction random_direction = ALL_CARDINALS[rng() % 4];
				command_queue.push_back(ship->move(random_direction));
			}
			else {
				command_queue.push_back(ship->stay_still());
			}
		}

		if (
			game.turn_number <= 200 &&
			me->halite >= constants::SHIP_COST &&
			!game_map->at(me->shipyard)->is_occupied())
		{
			command_queue.push_back(me->shipyard->spawn());
		}

		if (!game.end_turn(command_queue)) {
			break;
		}
	}

	return 0;
}

