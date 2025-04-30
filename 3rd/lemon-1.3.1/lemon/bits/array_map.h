/* -*- mode: C++; indent-tabs-mode: nil; -*-
 *
 * This file is a part of LEMON, a generic C++ optimization library.
 *
 * Copyright (C) 2003-2013
 * Egervary Jeno Kombinatorikus Optimalizalasi Kutatocsoport
 * (Egervary Research Group on Combinatorial Optimization, EGRES).
 *
 * Permission to use, modify and distribute this software is granted
 * provided that this copyright notice appears in all copies. For
 * precise terms see the accompanying LICENSE file.
 *
 * This software is provided "AS IS" with no warranty of any kind,
 * express or implied, and with no claim as to its suitability for any
 * purpose.
 *
 */

#ifndef LEMON_BITS_ARRAY_MAP_H
#define LEMON_BITS_ARRAY_MAP_H

#include <memory>

#include <lemon/bits/traits.h>
#include <lemon/bits/alteration_notifier.h>
#include <lemon/concept_check.h>
#include <lemon/concepts/maps.h>

 // \ingroup graphbits
 // \file
 // \brief Graph map based on the array storage.

namespace lemon {

    // \ingroup graphbits
    //
    // \brief Graph map based on the array storage.
    //
    // The ArrayMap template class is graph map structure that automatically
    // updates the map when a key is added to or erased from the graph.
    // This map uses the allocators to implement the container functionality.
    //
    // The template parameters are the Graph, the current Item type and
    // the Value type of the map.
    template <typename _Graph, typename _Item, typename _Value>
    class ArrayMap
        : public ItemSetTraits<_Graph, _Item>::ItemNotifier::ObserverBase {
    public:
        // The graph type.
        typedef _Graph GraphType;
        // The item type.
        typedef _Item Item;
        // The reference map tag.
        typedef True ReferenceMapTag;

        // The key type of the map.
        typedef _Item Key;
        // The value type of the map.
        typedef _Value Value;

        // The const reference type of the map.
        typedef const _Value& ConstReference;
        // The reference type of the map.
        typedef _Value& Reference;

        // The map type.
        typedef ArrayMap Map;

        // The notifier type.
        typedef typename ItemSetTraits<_Graph, _Item>::ItemNotifier Notifier;

    private:

        // The MapBase of the Map which implements the core registry function.
        typedef typename Notifier::ObserverBase Parent;

        typedef std::allocator<Value> Allocator;

    public:

        // \brief Graph initialized map constructor.
        explicit ArrayMap(const GraphType& graph) {
            Parent::attach(graph.notifier(Item()));
            allocate_memory();
            Notifier* nf = Parent::notifier();
            Item it;
            for (nf->first(it); it != INVALID; nf->next(it)) {
                int id = nf->id(it);
                std::allocator_traits<decltype(allocator)>::construct(allocator, &(values[id]), Value());
            }
        }

        // \brief Constructor to use default value to initialize the map.
        ArrayMap(const GraphType& graph, const Value& value) {
            Parent::attach(graph.notifier(Item()));
            allocate_memory();
            Notifier* nf = Parent::notifier();
            Item it;
            for (nf->first(it); it != INVALID; nf->next(it)) {
                int id = nf->id(it);
                std::allocator_traits<decltype(allocator)>::construct(allocator, &(values[id]), value);
            }
        }

    private:
        // \brief Constructor to copy a map of the same map type.
        ArrayMap(const ArrayMap& copy) : Parent() {
            if (copy.attached()) {
                attach(*copy.notifier());
            }
            capacity = copy.capacity;
            if (capacity == 0) return;
            values = allocator.allocate(capacity);
            Notifier* nf = Parent::notifier();
            Item it;
            for (nf->first(it); it != INVALID; nf->next(it)) {
                int id = nf->id(it);
                std::allocator_traits<decltype(allocator)>::construct(allocator, &(values[id]), copy.values[id]);
            }
        }

        // \brief Assign operator.
        ArrayMap& operator=(const ArrayMap& cmap) {
            return operator=<ArrayMap>(cmap);
        }

        // \brief Template assign operator.
        template <typename CMap>
        ArrayMap& operator=(const CMap& cmap) {
            checkConcept<concepts::ReadMap<Key, _Value>, CMap>();
            const typename Parent::Notifier* nf = Parent::notifier();
            Item it;
            for (nf->first(it); it != INVALID; nf->next(it)) {
                set(it, cmap[it]);
            }
            return *this;
        }

    public:
        // \brief The destructor of the map.
        virtual ~ArrayMap() {
            if (attached()) {
                clear();
                detach();
            }
        }

    protected:
        using Parent::attach;
        using Parent::detach;
        using Parent::attached;

    public:

        // \brief The subscript operator.
        Value& operator[](const Key& key) {
            int id = Parent::notifier()->id(key);
            return values[id];
        }

        // \brief The const subscript operator.
        const Value& operator[](const Key& key) const {
            int id = Parent::notifier()->id(key);
            return values[id];
        }

        // \brief Setter function of the map.
        void set(const Key& key, const Value& val) {
            (*this)[key] = val;
        }

    protected:

        // \brief Adds a new key to the map.
        virtual void add(const Key& key) {
            Notifier* nf = Parent::notifier();
            int id = nf->id(key);
            if (id >= capacity) {
                int new_capacity = (capacity == 0 ? 1 : capacity);
                while (new_capacity <= id) {
                    new_capacity <<= 1;
                }
                Value* new_values = allocator.allocate(new_capacity);
                Item it;
                for (nf->first(it); it != INVALID; nf->next(it)) {
                    int jd = nf->id(it);
                    if (id != jd) {
                        std::allocator_traits<decltype(allocator)>::construct(allocator, &(new_values[jd]), values[jd]);
                        std::allocator_traits<decltype(allocator)>::destroy(allocator, &(values[jd]));
                    }
                }
                if (capacity != 0) allocator.deallocate(values, capacity);
                values = new_values;
                capacity = new_capacity;
            }
            std::allocator_traits<decltype(allocator)>::construct(allocator, &(values[id]), Value());
        }

        // \brief Adds more new keys to the map.
        virtual void add(const std::vector<Key>& keys) {
            Notifier* nf = Parent::notifier();
            int max_id = -1;
            for (int i = 0; i < int(keys.size()); ++i) {
                int id = nf->id(keys[i]);
                if (id > max_id) {
                    max_id = id;
                }
            }
            if (max_id >= capacity) {
                int new_capacity = (capacity == 0 ? 1 : capacity);
                while (new_capacity <= max_id) {
                    new_capacity <<= 1;
                }
                Value* new_values = allocator.allocate(new_capacity);
                Item it;
                for (nf->first(it); it != INVALID; nf->next(it)) {
                    int id = nf->id(it);
                    bool found = false;
                    for (int i = 0; i < int(keys.size()); ++i) {
                        int jd = nf->id(keys[i]);
                        if (id == jd) {
                            found = true;
                            break;
                        }
                    }
                    if (found) continue;
                    std::allocator_traits<decltype(allocator)>::construct(allocator, &(new_values[id]), values[id]);
                    std::allocator_traits<decltype(allocator)>::destroy(allocator, &(values[id]));
                }
                if (capacity != 0) allocator.deallocate(values, capacity);
                values = new_values;
                capacity = new_capacity;
            }
            for (int i = 0; i < int(keys.size()); ++i) {
                int id = nf->id(keys[i]);
                std::allocator_traits<decltype(allocator)>::construct(allocator, &(values[id]), Value());
            }
        }

        // \brief Erase a key from the map.
        virtual void erase(const Key& key) {
            int id = Parent::notifier()->id(key);
            std::allocator_traits<decltype(allocator)>::destroy(allocator, &(values[id]));
        }

        // \brief Erase more keys from the map.
        virtual void erase(const std::vector<Key>& keys) {
            for (int i = 0; i < int(keys.size()); ++i) {
                int id = Parent::notifier()->id(keys[i]);
                std::allocator_traits<decltype(allocator)>::destroy(allocator, &(values[id]));
            }
        }

        // \brief Builds the map.
        virtual void build() {
            Notifier* nf = Parent::notifier();
            allocate_memory();
            Item it;
            for (nf->first(it); it != INVALID; nf->next(it)) {
                int id = nf->id(it);
                std::allocator_traits<decltype(allocator)>::construct(allocator, &(values[id]), Value());
            }
        }

        // \brief Clear the map.
        virtual void clear() {
            Notifier* nf = Parent::notifier();
            if (capacity != 0) {
                Item it;
                for (nf->first(it); it != INVALID; nf->next(it)) {
                    int id = nf->id(it);
                    std::allocator_traits<decltype(allocator)>::destroy(allocator, &(values[id]));
                }
                allocator.deallocate(values, capacity);
                capacity = 0;
            }
        }

    private:

        void allocate_memory() {
            Notifier* nf = Parent::notifier();
            Item it;
            int max_id = -1;
            for (nf->first(it); it != INVALID; nf->next(it)) {
                int id = nf->id(it);
                if (id > max_id) {
                    max_id = id;
                }
            }
            capacity = max_id + 1;
            values = allocator.allocate(capacity);
        }

        // The allocator used for storing the values.
        Allocator allocator;

        // The values in the map.
        Value* values;

        // The capacity of the map.
        int capacity;
    };

} // namespace lemon

#endif /* LEMON_BITS_ARRAY_MAP_H */
