// CPCS 324 Algorithms & Data Structures 2
// Outline - Priority queue data structure
// 2019, Dr. Muhammad Al-Hashimi


// -----------------------------------------------------------------------
// Basic design decisions and implementation planning (objects & interfaces)

// initial requirements: to quickly support Dijkstra and second Prim's algorithms, 
// implement minimum priority functionality

// design decisions:
// Reuse the 324 linked list implementation
// how will the PQ ops be implemented?
// <student fill here>

// code plan: start to write your API specifications (JSDOC comments) and think about 
// design consequences (design impact)

// Impact analysis:
// <student fill here>
// -----------------------------------------------------------------------

// Priority queue object constructor (document using JSDOC comments)

/**
 * @description Priority queue object constructor
 * @method PQueue
 * @author Esraa ALZahrani
 */
function PQueue() {
    this.pq = new List(); // requirement: linked-list implementation
    // specify (design) methods
    this.isEmpty = isEmptyImp; // return true if queue empty
    this.deleteMin = deleteMinImp; // remove/return item with minimum priority
    this.insert = insertImp; // insert/update an item with priority

}

/**
 * @description Priority queue node constructor 
 * @method PQnode
 * @param {number} item
 * @param {number} key
* @author Esraa ALZahrani
 */
// -----------------------------------------------------------------------
// Priority queue node constructor (document using JSDOC comments)
function PQNode(item, key) {
    this.item = item;
    this.prior = key;

    // specify (design) methods

}
// -----------------------------------------------------------------------
// functions used by PQueue() object methods
// specify interface information (JSDOC comments)
// function names should not clash with linklist.js and queue.js



//~~~~~~~~~~~
/** 
 * @method isEmpty 
* @description check if PQ is empty
* @returns {boolean}
* @author Esraa ALZahrani
* */
function isEmptyImp() {
    return (this.pq.isEmpty());
}


//~~~~~~~~~~~
/**
 * @description deletes PQ node with highest priority (smaller number = higher priority)
 * @method deleteMin
 * @author Esraa ALZahrani
 * @returns {PQNode} node deleted
 */
function deleteMinImp() {
    if (!isEmpty()) {
        var Ptr = this.pq.first;
        var temp = ptr.next;
        while (Ptr != null) {
            if (temp.prior > Ptr.item.prior) {
                temp = Ptr;
                break;
            }
            ptr = Ptr.next;
        }
        var cur = temp.item;
        this.temp = this.temp.next;
        return cur;
    }
}

//~~~~~~~~~~~
/**
* 
* @description implement function to add new node or update the priority of the item is oready exist
* @param {number} item 
* @param {number} key 
*/
function insertImp(item, key) {
    var pqnode = new PQNode(item, key);
    var Ptr = this.pq.first;
    var exist;
    if (this.isEmpty()) {
        this.pq.insert(pqnode);
    }
    else {
        while (Ptr != null) {
            if (pqnode.item === Ptr.item) {
                Ptr.item.prior = key;
                exist = true;
                break;
            }
            ptr = Ptr.next;
        }
        if (!exist) {
            this.pq.insert(pqnode);
        }
    }
}